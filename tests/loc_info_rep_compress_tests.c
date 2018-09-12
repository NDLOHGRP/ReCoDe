/*
Tests the compressibility of location information: 2-D spatial coordinates of pixel accurate puddle centroids,
using four data representations: 1. 24-bit xy coordinates, 2. 24-bit xy coordinates sorted by x, 3. 24-bit xy coordinates represented as quad-tree and 4. binary image representation
and two compression algorithms: 1. g-zip optimized for speed and 2. g-zip optimized for compression
*/

/*
Representation 1: 24-bit xy coordinates
	1. Read data image and dark image
	2. Perform dark subtraction (thresholding)
	3. Perform connected components analysis
	4. Compute centroids
	5. Pack 12-bit pairs into 3 bytes
	
Representation 2: 24-bit xy coordinates sorted by X
	1. Read data image and dark image
	2. Perform dark subtraction (thresholding)
	3. Perform connected components analysis
	4. Compute centroids
	5. Sort centroids by X
	6. Pack 12-bit pairs into 3 bytes
	
Representation 3: Quad-tree
	1. Read data image and dark image
	2. Perform dark subtraction (thresholding)
	3. Perform connected components analysis
	4. Compute centroids
	5. For each centroid find: binary-tree it belongs to and offset within binary-tree
	6. Pack binary tree counts
	7. Pack offsets

Representation 4: Binary image
	1. Read data image and dark image
	2. Perform dark subtraction (thresholding)
	3. Perform connected components analysis
	4. Compute centroids
	5. Make binary image by dressing centroids
*/




/*****************************************************************************************
to compile:	gcc loc_info_rep_compress_tests.c -lz -O3 -o loc_info_rep_compress_tests -lm
to run:		./loc_info_rep_compress_tests

-lz: links zlib
-lm: links math.h NOTE: it is important to keep -lm as the last option
-o:	 write build output to file
-O3: optimizer flag
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "zlib.h"

#define X_LEN 4096
#define Y_LEN 4096
#define MAX_ELEM (X_LEN*Y_LEN) 
#define MAX_CLUSTER (X_LEN*Y_LEN) 
#define WORDLEN 300

#define  NUMFRAMES		1
#define  FRAMEXDIM  	4096
#define  FRAMEYDIM  	4096
#define  FRAMESIZE		FRAMEXDIM * FRAMEYDIM
#define	 DTYPE_USHORT	16

#define SetBit(A,k)		( A[(k/8)] |= (1 << (k%8)) )
#define ClearBit(A,k)	( A[(k/8)] &= ~(1 << (k%8)) )
#define CheckBit(A,k)	( A[(k/8)] & (1 << (k%8)) )

typedef struct {
	unsigned short	nx;
	unsigned short	ny;
	unsigned long	nz;
	unsigned char	dtype;
} DataSize;


int GetMaxCompressedLen(int nLenSrc) {
	int n16kBlocks = (nLenSrc + 16383) / 16384; // round up any fraction of a block
	return (nLenSrc + 6 + (n16kBlocks * 5));
}

char* concat(const char *s1, const char *s2) {
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = (char*) malloc(len1+len2+1);//+1 for the null-terminator
    //check for errors in malloc
    memcpy(result, s1, len1);
    memcpy(result+len1, s2, len2+1);//+1 to copy the null-terminator
    return result;
}

void writeToFile (const char* filename, unsigned char* data, int dataLength) {
	
	FILE * fp;
	fp = fopen (filename, "w+");
	
	int i = 0;
	for (i = 0; i < dataLength; i++) {
		fprintf(fp, "%u ", data[i]);
	}
	
	fclose(fp);
	
}

void saveFrame (const char* filename, unsigned short* buffer, DataSize d, int framenum) {
	
	int r = 1;
	int c = 1;
	
	FILE * fp;

	fp = fopen (filename, "w+");
	
	long int offset = d.ny*d.nx*framenum;
	for (r = 0; r < d.ny; r++) {
		for (c = 0; c < d.nx; c++) {
			fprintf(fp, "%d ", buffer[offset + r*d.nx + c]);
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
}

float gzip_compress_stream (unsigned char* data, unsigned long nDataLen, int mode, long *nCompressedSize, unsigned char** compressedData) {
	
	clock_t p_start = clock();
	
	unsigned int nLenDst = GetMaxCompressedLen(nDataLen);

	//unsigned char* pbDst = new unsigned char[nLenDst];  // alloc dest buffer
	unsigned char* pbDst = (unsigned char*)calloc(nLenDst, 1);
	
	// zlib struct
    z_stream zInfo;
    zInfo.zalloc	= Z_NULL;
    zInfo.zfree 	= Z_NULL;
    zInfo.opaque 	= Z_NULL;
	
    // setup input and compressed output
    zInfo.avail_in 	= nDataLen; 		// size of input
	zInfo.total_in	= nDataLen;
    zInfo.avail_out = nLenDst; 			// size of output
	zInfo.total_out	= nLenDst;
	zInfo.next_in 	= (Bytef *)data; 	// input char array
    zInfo.next_out 	= (Bytef *)pbDst; 	// output char array
	
	// the actual compression work
    int nErr = -1;
	*nCompressedSize = -1;
	if (mode == 0) {
		nErr = deflateInit(&zInfo, Z_BEST_SPEED); 		// zlib function
	} else if(mode == 1) {
		nErr = deflateInit(&zInfo, Z_BEST_COMPRESSION); // zlib function
	}
	if (nErr == Z_OK) {
		nErr = deflate(&zInfo, Z_FINISH);              	// zlib function
		if (nErr == Z_STREAM_END) {
			*nCompressedSize = zInfo.total_out;
		}
	}
	deflateEnd(&zInfo);    // zlib function
	
	*compressedData = pbDst;
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
	
}

float gzip_decompress_stream (unsigned char* compressedData, unsigned char* deCompressedData, unsigned int nCompressedLen, unsigned long nDataLen) {
	
	clock_t p_start = clock();
	
	// zlib struct
    z_stream zInfo;
    zInfo.zalloc	= Z_NULL;
    zInfo.zfree 	= Z_NULL;
    zInfo.opaque 	= Z_NULL;
	
	zInfo.total_in	= nCompressedLen;
	zInfo.avail_in	= nCompressedLen;
    zInfo.total_out	= nDataLen;
	zInfo.avail_out	= nDataLen;
    zInfo.next_in	= (Bytef *)compressedData;
    zInfo.next_out	= deCompressedData;

	printf("%d\n", nCompressedLen);
	printf("%d\n", nDataLen);
	printf("%u\n", compressedData[0]);
	
    int nErr, nRet= -1;
    nErr = inflateInit( &zInfo );            // zlib function
    if (nErr == Z_OK) {
		printf("In 1.\n");
        nErr = inflate(&zInfo, Z_FINISH);   // zlib function
        if (nErr == Z_STREAM_END) {
			printf("In 2.\n");
            nRet = zInfo.total_out;
        }
    }
    inflateEnd (&zInfo);   					// zlib function
    printf("De-compressed Size: %d\n", nRet);
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
	
}

void serialize_L4 (const char* filename, unsigned char* compressedData, int nCompressedSize, DataSize d) {
	
	unsigned char level = 4;
	unsigned char decimation = 0;
	unsigned short nx = d.nx;
	unsigned short ny = d.ny;
	unsigned short nz = d.nz;
	
	FILE *fp = fopen (filename, "wb");
	fwrite (&level, sizeof(char), 1, fp);
	fwrite (&decimation, sizeof(char), 1, fp);
	fwrite (&nx, sizeof(short), 1, fp);
	fwrite (&ny, sizeof(short), 1, fp);
	fwrite (&nz, sizeof(short), 1, fp);
	fwrite (compressedData , sizeof(char), nCompressedSize, fp);
	fclose (fp);
}

void load_L4 (const char* filename, unsigned char* compressedData, int* nCompressedSize, DataSize* d) {
	
	
}

float threshold_image (unsigned short* frameBuffer, unsigned short* darkBuffer, unsigned short epsilon, DataSize h, int *clusterMap) {
	
	clock_t p_start = clock();
	
	unsigned int row, col;
	unsigned long linear_index;
	unsigned short tmp_f, tmp_d, thresh_xy;
	
	//Label the clusters as -1, nonclusters 0
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			tmp_d = darkBuffer[linear_index];
			thresh_xy = tmp_d + epsilon;
			clusterMap[linear_index] = (tmp_f > thresh_xy) ? -1 : 0;
		}
	}

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}

// Used in L3
float threshold_and_make_binary_image (unsigned short* frameBuffer, unsigned short* darkBuffer, unsigned short epsilon, DataSize h, unsigned char *binaryImage) {
	
	clock_t p_start = clock();
	
	unsigned int row, col;
	unsigned long linear_index;
	unsigned short tmp_f, tmp_d, thresh_xy;
	
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			tmp_d = darkBuffer[linear_index];
			thresh_xy = tmp_d + epsilon;

			if (tmp_f > thresh_xy) {
				SetBit(binaryImage, linear_index);
			}
		}
	}

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}

// Used in L1
float threshold_make_binary_image_and_get_pixvals (unsigned short* frameBuffer, unsigned short* darkBuffer, unsigned short epsilon, DataSize h, unsigned char *binaryImage, unsigned short *pixvals, unsigned long* nFgPixels, unsigned short *data_min, unsigned short* data_max) {
	
	clock_t p_start = clock();
	
	unsigned int row, col;
	unsigned long linear_index;
	unsigned short tmp_f, tmp_d, thresh_xy, dark_subtracted_pixval;
	
	*nFgPixels = 0;
	*data_min = 65535;
	*data_max = 0;
	
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			tmp_d = darkBuffer[linear_index];
			thresh_xy = tmp_d + epsilon;

			if (tmp_f > thresh_xy) {
				SetBit(binaryImage, linear_index);
				dark_subtracted_pixval = tmp_f - thresh_xy;
				pixvals[(*nFgPixels)++] = dark_subtracted_pixval;
				if (dark_subtracted_pixval > *data_max) {
					*data_max = dark_subtracted_pixval;
					printf("MAX: %hu\n", *data_max);
				}
				if (dark_subtracted_pixval < *data_min) {
					*data_min = dark_subtracted_pixval;
				}
			} else {
				*data_min = 0;
			}
		}
	}

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}

// Used for validation
float get_dark_subtracted_image (unsigned short* frameBuffer, unsigned short* darkBuffer, unsigned short epsilon, DataSize h, unsigned short *darkSubtractedImage) {
	
	clock_t p_start = clock();
	
	unsigned int row, col;
	unsigned long linear_index;
	unsigned short tmp_f, tmp_d, thresh_xy;
	
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			tmp_d = darkBuffer[linear_index];
			thresh_xy = tmp_d + epsilon;

			if (tmp_f > thresh_xy) {
				darkSubtractedImage[linear_index] = tmp_f - thresh_xy;
			} else {
				darkSubtractedImage[linear_index] = 0;
			}
		}
	}

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}

float forest_fire_conn_comp_labelling (unsigned short* frameBuffer, unsigned short* darkBuffer, unsigned short epsilon, DataSize h, int *clusterMap, float *x_coor, float *y_coor, int *label) {
	
	clock_t p_start = clock();
	
	unsigned int x, y, x0, y0, x1, y1;
	unsigned long l0, l1, l2, l_curr;
	
	//set values of cluster map boundary to zero
	for (x = 0; x < h.nx; x++) {
		y0 = 0;
		y1 = Y_LEN - 1;
		l0 = y0 * h.nx + x;
		l1 = y1 * h.nx + x;
		clusterMap[l0] = 0;
		clusterMap[l1] = 0;
	}
	
	for (y = 0; y < h.ny; y++) {
		x0 = 0;
		x1 = h.nx - 1;
		l0 = y * h.nx + x0;
		l1 = y * h.nx + x1;
		clusterMap[l0] = 0;
		clusterMap[l1] = 0;
	}
	
	int n_max = 4;
	int nb[4][2] = { { -1,0 },{ 1,0 },{ 0,-1 },{ 0,1 } };

	//Identify and label clusters
	int *bd_x, *bd_y;					// search list for processing current cluster
	bd_x = (int *) malloc(MAX_ELEM * sizeof(*bd_x));
	bd_y = (int *) malloc(MAX_ELEM * sizeof(*bd_y));
	
	unsigned short dark_subtracted_pixel_val;
	
	for (y = 0; y < h.ny; y++)
	{
		for (x = 0; x < h.nx; x++)
		{
			int l1, l2, l_curr;
			l1 = y * h.nx + x;
			if (clusterMap[l1] == -1)
			{
				int n;
				int n_elem = 0; 	// pointer to top of search list
				int weight_sum = 0;
				unsigned int x_curr, y_curr, x2, y2;

				// init search list for a new cluster`
				bd_x[n_elem] = x;
				bd_y[n_elem] = y;
				
				// set new cluster label
				(*label)++;
				
				do {
					//process the next pixel in search list
					x_curr = bd_x[n_elem];
					y_curr = bd_y[n_elem];
					l_curr = y_curr * h.nx + x_curr;

					clusterMap[l_curr] = *label;

					dark_subtracted_pixel_val = frameBuffer[l_curr] - (darkBuffer[l_curr] + epsilon);
					dark_subtracted_pixel_val = (dark_subtracted_pixel_val > 0) ? dark_subtracted_pixel_val : 0;
					
					x_coor[*label] += x_curr*dark_subtracted_pixel_val;
					y_coor[*label] += y_curr*dark_subtracted_pixel_val;
					weight_sum += dark_subtracted_pixel_val;
				
					//printf("%.2f, %.2f\n", x_coor[*label], y_coor[*label]);
				
					n_elem--;	// pop one

					//check current pixel's neighbors and add them to search list if they are fg
					for (n = 0; n < n_max; n++)
					{
						x2 = x_curr + nb[n][0];
						y2 = y_curr + nb[n][1];
						l2 = y2 * h.nx + x2;
						if (clusterMap[l2] == -1)
						{
							n_elem++;
							bd_x[n_elem] = x2;
							bd_y[n_elem] = y2;
							clusterMap[l2] = *label;
						}
					}

				} while (n_elem >= 0);	// continue processing this cluster till search list is empty

				//printf("Processed Clusters = %d\n", *label);
				x_coor[*label] /= (weight_sum*1.0);
				y_coor[*label] /= (weight_sum*1.0);
				//printf("%.2f, %.2f, %d\n", x_coor[*label], y_coor[*label], weight_sum);
				//printf("%d, %d\n", x, y);
			}
		}
	}
	
	free(bd_x);
	free(bd_y);
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
	
}

float scale_and_pack_pixvals(unsigned short *pixvals, unsigned long nFgPixels, unsigned short data_min, unsigned short data_max, unsigned char bit_depth, unsigned char *packedPixvals) {
	
	int p, n, linear_index, nth_bit_of_pth_pixval;
	unsigned short scaled_pixval;
	double pixval_01;
	
	unsigned short MAX_VAL = pow(2, bit_depth) - 1;
	
	for (p=0; p<nFgPixels; p++) {
		for (n=0; n<bit_depth; n++) {
			
			pixval_01 = ((pixvals[p] - data_min)*1.0)/((data_max - data_min)*1.0);
			scaled_pixval = (unsigned short)ceil(pixval_01*MAX_VAL);
			
			nth_bit_of_pth_pixval = (scaled_pixval & ( 1 << n )) >> n;
			if (nth_bit_of_pth_pixval == 1) {
				linear_index = p*bit_depth + n;
				SetBit(packedPixvals, linear_index);
			}
			//printf("%f, %hu, %hu, %hu, %hu\n", pixval_01, scaled_pixval, pixvals[p], data_min, data_max);
		}
	}
	
	printf("MAX_VAL for given bit_depth: %hu\n", MAX_VAL);
	
}

void do_L1_Reduce_Compress (unsigned short* frameBuffer, unsigned short* darkBuffer, int epsilon, unsigned long n_Frames, int mode, unsigned char bit_depth, float* compress_times_sizes) {

	float thresh_time, packing_time, reduction_time, compression_time, compression_time_1, compression_time_2;

	printf("Size of a short is: %d\n", sizeof(short));
	printf("Size of a char is: %d\n", sizeof(char));

	DataSize h = { 4096, 4096, n_Frames, 16 };
	unsigned char *binaryImage, *packedPixvals;
	unsigned short* pixvals;
	char epsilon_s = (unsigned short)epsilon;
	
	printf("There are %d frames. Starting to cluster on them..\n", h.nz);
	
	unsigned char part_id = 1;
	char* part_num = malloc(sizeof(char)*(int)log10(part_id));
	sprintf(part_num, "%d", part_id);
	const char* filename = concat(concat("Part_", part_num), ".rc1");
	free(part_num);
	
	unsigned char reduction_level = 1;
	FILE *fp = fopen (filename, "wb");
	fwrite (&part_id, sizeof(char), 1, fp);
	fwrite (&reduction_level, sizeof(char), 1, fp);
	fwrite (&bit_depth, sizeof(char), 1, fp);
	fwrite (&h.nx, sizeof(short), 1, fp);
	fwrite (&h.ny, sizeof(short), 1, fp);
	fwrite (&h.nz, sizeof(long), 1, fp);
	
	unsigned int z;
	for (z = 0; z < h.nz; z++) {
		
		unsigned long nFgPixels = 0;
		unsigned short data_min = 65535;
		unsigned short data_max = 0;
		
		pixvals = (unsigned short*)calloc((h.nx * h.ny), sizeof(short));			// allocated max space - where all pixels are fg pixels, must calloc to init to 0
		binaryImage = (unsigned char*)calloc((h.nx * h.ny) / 8, sizeof(char));		// unlike clusterMap for L4, binaryImage has to be re-inited for every frame here
		
		thresh_time = threshold_make_binary_image_and_get_pixvals (frameBuffer, darkBuffer, epsilon_s, h, binaryImage, pixvals, &nFgPixels, &data_min, &data_max);
		
		unsigned long bytesRequiredForPacking = ceil((nFgPixels*bit_depth)/8);
		packedPixvals = (unsigned char*)calloc(bytesRequiredForPacking, sizeof(char));
		
		packing_time = scale_and_pack_pixvals(pixvals, nFgPixels, data_min, data_max, bit_depth, packedPixvals);
		
		reduction_time = thresh_time + packing_time;
		
		unsigned long nLenOrig_BinaryImage = (h.nx * h.ny) / 8;
		unsigned long nCompressedSize, nCompressedSize_1, nCompressedSize_2;
		unsigned char* compressedBinaryImage;
		unsigned char* compressedPixvals;
		compression_time_1 = gzip_compress_stream (binaryImage, nLenOrig_BinaryImage, mode, &nCompressedSize_1, &compressedBinaryImage);
		compression_time_2 = gzip_compress_stream (packedPixvals, bytesRequiredForPacking, mode, &nCompressedSize_2, &compressedPixvals);
		
		compression_time = compression_time_1 + compression_time_2;
		nCompressedSize = nCompressedSize_1 + nCompressedSize_2;
		
		printf("Compression done.\n");
		//printf("Max Compressed Size: %d\n", nLenDst);
		printf("Compressed Size: %d\n", nCompressedSize);
		printf("Compression time for frame %d: %f ms\n", z, compression_time);

		printf("Bytes Required for packing: %d\n", bytesRequiredForPacking);
		fwrite (&nCompressedSize_1, sizeof(long), 1, fp);
		fwrite (&nCompressedSize_2, sizeof(long), 1, fp);
		fwrite (&bytesRequiredForPacking, sizeof(long), 1, fp);
		fwrite (compressedBinaryImage , sizeof(char), nCompressedSize_1, fp);
		fwrite (compressedPixvals , sizeof(char), nCompressedSize_2, fp);
		
		compress_times_sizes[0] = reduction_time;
		compress_times_sizes[1] = compression_time;
		compress_times_sizes[2] = (h.nx * h.ny)*1.5;
		compress_times_sizes[3] = nCompressedSize*1.0;
		
		/*
		unsigned char* deCompressedBinaryImage = (unsigned char*)calloc(nLenOrig_BinaryImage, sizeof(char));
		unsigned char* deCompressedPixvals = (unsigned char*)calloc(bytesRequiredForPacking, sizeof(char));
		
		gzip_decompress_stream(compressedBinaryImage, deCompressedBinaryImage, nCompressedSize_1, nLenOrig_BinaryImage);
		gzip_decompress_stream(compressedPixvals, deCompressedPixvals, nCompressedSize_2, bytesRequiredForPacking);
		
		writeToFile("original_binary_image.txt", binaryImage, nLenOrig_BinaryImage);
		writeToFile("de-compressed_binary_image.txt", deCompressedBinaryImage, nLenOrig_BinaryImage);
		
		printf("Original Data: %u, %u, %u\n", binaryImage[0], binaryImage[10], binaryImage[20]);
		printf("De-compressed Data: %u, %u, %u\n", deCompressedBinaryImage[0], deCompressedBinaryImage[10], deCompressedBinaryImage[20]);
		*/
		
		writeToFile("original_binary_image.txt", binaryImage, nLenOrig_BinaryImage);
		
		DataSize d = { 4096, 4096, n_Frames, 16 };
		unsigned short *darkSubtractedImage = (unsigned short*)calloc((h.nx * h.ny), sizeof(short));
		get_dark_subtracted_image (frameBuffer, darkBuffer, epsilon_s, d, darkSubtractedImage);
		saveFrame ("original_dark_subtracted_image.txt", darkSubtractedImage, d, 0);
		
		free(pixvals);
		free(packedPixvals);
		free(binaryImage);
		free(compressedBinaryImage);
		free(compressedPixvals);
		
		//free(deCompressedBinaryImage);
		//free(deCompressedPixvals);
	}

	fclose (fp);
	
}

void do_L1_Decompress_Expand (const char* compressed_filename) {
	
	//, unsigned short** frameBuffer
	
	FILE* fp = fopen(compressed_filename,"rb");
	
	unsigned char part_id;
	fread (&part_id, sizeof(char), 1, fp);
	printf("Part ID: %u\n", part_id);
	
	unsigned char reduction_level;
	fread (&reduction_level, sizeof(char), 1, fp);
	printf("Reduction Level: %u\n", reduction_level);
	
	unsigned char bit_depth;
	fread (&bit_depth, sizeof(char), 1, fp);
	printf("Bit Depth: %u\n", bit_depth);
	
	unsigned short n_cols;
	fread (&n_cols, sizeof(short), 1, fp);
	printf("Num Columns: %u\n", n_cols);
	
	unsigned short n_rows;
	fread (&n_rows, sizeof(short), 1, fp);
	printf("Num Rows: %u\n", n_rows);
	
	unsigned long n_Frames;
	fread (&n_Frames, sizeof(long), 1, fp);
	printf("Num Frames: %d\n", n_Frames);
	
	unsigned int sz = n_cols*n_rows*n_Frames;
	printf("sz = %d\n", sz);
	unsigned short* frameBuffer = (unsigned short*)calloc(sz, sizeof(short));
	printf("JFC: %u\n", frameBuffer[0]);
	
	int z;
	for (z = 0; z < 1; z++) {
		
		unsigned long nCompressedSize_BinaryImage;
		fread (&nCompressedSize_BinaryImage, sizeof(long), 1, fp);
		printf("Compressed Size 1: %d\n", nCompressedSize_BinaryImage);
		
		unsigned long nCompressedSize_Pixvals;
		fread (&nCompressedSize_Pixvals, sizeof(long), 1, fp);
		printf("Compressed Size 2: %d\n", nCompressedSize_Pixvals);
		
		unsigned long bytesRequiredForPacking;
		fread (&bytesRequiredForPacking, sizeof(long), 1, fp);
		printf("Packed Bytes: %d\n", bytesRequiredForPacking);
		
		unsigned char* compressedBinaryImage = (unsigned char*)calloc(nCompressedSize_BinaryImage, sizeof(char));
		fread (compressedBinaryImage, sizeof(char), nCompressedSize_BinaryImage, fp);
		
		unsigned char* compressedPixvals = (unsigned char*)calloc(nCompressedSize_Pixvals, sizeof(char));
		fread (compressedPixvals, sizeof(char), nCompressedSize_Pixvals, fp);
		
		unsigned long nLenOrig_BinaryImage = (n_cols * n_rows) / 8;
		unsigned char* deCompressedBinaryImage = (unsigned char*)calloc(nLenOrig_BinaryImage, sizeof(char));
		unsigned char* deCompressedPixvals = (unsigned char*)calloc(bytesRequiredForPacking, sizeof(char));
		
		gzip_decompress_stream(compressedBinaryImage, deCompressedBinaryImage, nCompressedSize_BinaryImage, nLenOrig_BinaryImage);
		gzip_decompress_stream(compressedPixvals, deCompressedPixvals, nCompressedSize_Pixvals, bytesRequiredForPacking);

		writeToFile("de-compressed_binary_image.txt", deCompressedBinaryImage, nLenOrig_BinaryImage);
		
		printf("bytesRequiredForPacking: %d\n", bytesRequiredForPacking);
		
		unsigned long count = 0;
		unsigned long frame_offset = n_cols * n_rows * z;
		unsigned long bit_offset = 0;
		unsigned long pixel_offset = 0;
		unsigned short pixval = 0;

		int r,c;
		for (r=0; r<n_rows; r++) {
			for (c=0; c<n_cols; c++) {				
				pixel_offset = r*n_cols + c;
				//printf("%d, %d, %lu\n", r, c, pixel_offset);
				if (CheckBit(deCompressedBinaryImage, pixel_offset)) {
					bit_offset = bit_depth*count;
					//printf("%lu, %lu, %lu\n", count, pixel_offset, bit_offset);

					pixval = 0;
					unsigned int n;
					for (n=0; n<bit_depth; n++) {
						unsigned long off = bit_offset + n;
						if (CheckBit(deCompressedPixvals, off)) {
							pixval += pow(2,n);
						}
					}

					frameBuffer[frame_offset + pixel_offset] = pixval;
					count++;

				}
			}
		}

		DataSize d = { n_cols, n_rows, n_Frames, 16 };
		saveFrame ("de-compressed_original_image.txt", frameBuffer, d, 0);
		
		free(compressedBinaryImage);
		free(compressedPixvals);
		free(deCompressedBinaryImage);
		free(deCompressedPixvals);
		free(frameBuffer);
		
		printf("Here.\n");
		
	}
	
	fclose(fp);
}

void do_L4_Decompress_Expand_Binary_Image (unsigned char* compressedBuffer, unsigned short* outBuffer, DataSize d, float* decompress_times_sizes) {
	
}

void do_L4_Reduce_Compress_Binary_Image (unsigned short* frameBuffer, unsigned short* darkBuffer, int epsilon, int n_Frames, float* compress_times_sizes, int mode) {

	float thresh_time, cca_time, reduction_time, compression_time;	

	int fsize = sizeof(float), isize = sizeof(int), csize = sizeof(char), lsize = sizeof(long);

	//printf("Size of a int is: %d\n", lsize);

	DataSize h = { 4096, 4096, n_Frames, 16 };
	int* clusterMap = (int *) malloc(h.nx * h.ny * isize);
	unsigned char *centroidImage;
	char epsilon_s = (unsigned short)epsilon;
	
	printf("There are %d frames. Starting to cluster on them..\n", h.nz);

	unsigned int z;
	for (z = 0; z < h.nz; z++) {
		
		thresh_time = threshold_image (frameBuffer, darkBuffer, epsilon_s, h, clusterMap);
		
		printf("Background subtraction done.\n");
		
		float *x_coor, *y_coor;
		int label = 0;
		x_coor = (float *) malloc(MAX_CLUSTER * sizeof(*x_coor));
		y_coor = (float *) malloc(MAX_CLUSTER * sizeof(*y_coor));
		cca_time = forest_fire_conn_comp_labelling (frameBuffer, darkBuffer, epsilon_s, h, clusterMap, x_coor, y_coor, &label);
		
		
		unsigned int x;
		unsigned long l;
		centroidImage = (unsigned char*)calloc((h.nx * h.ny) / 8, 1);
		for (x = 1; x <= label; x++) {
			l = (int)round(y_coor[x]) * h.nx + (int)round(x_coor[x]);
			SetBit(centroidImage, l);
		}
		printf("No. of Labels: %d\n", label);
		reduction_time = thresh_time + cca_time;
		
		unsigned long nLenOrig = (h.nx * h.ny) / 8;
		unsigned long nCompressedSize;
		unsigned char* compressedCentroidImage;
		compression_time = gzip_compress_stream (centroidImage, nLenOrig, mode, &nCompressedSize, &compressedCentroidImage);
		
		printf("Compression done.\n");
		//printf("Max Compressed Size: %d\n", nLenDst);
		printf("Compressed Size: %d\n", nCompressedSize);
		printf("Compression time for frame %d: %f ms\n", z, compression_time);

		compress_times_sizes[0] = reduction_time;
		compress_times_sizes[1] = compression_time;
		compress_times_sizes[2] = nLenOrig*1.0;
		compress_times_sizes[3] = nCompressedSize*1.0;
		
		free(x_coor);
		free(y_coor);
		//free(max);
	}

	free(clusterMap);
	free(centroidImage);
	
}

void do_L4_Reduce_Compress_Binary_Image_Old (unsigned short* frameBuffer, unsigned short* darkBuffer, int epsilon, int n_Frames, float* compress_times_sizes, int mode) {

	char epsilon_s = (unsigned short)epsilon;

	int error;
	int x, y, z;
	unsigned int l;

	unsigned short dark_subtracted_pixel_val;

	unsigned char *centroidImage;

	int *clusterMap;
	int tmp_f, tmp_d, thresh_xy;

	clock_t reduction_start, reduction_end, compression_start, compression_end;
	float reduction_time, compression_time;

	float *x_coor, *y_coor;
	int *bd_x, *bd_y;

	x_coor = (float *) malloc(MAX_CLUSTER * sizeof(*x_coor));
	y_coor = (float *) malloc(MAX_CLUSTER * sizeof(*y_coor));

	int fsize = sizeof(float), isize = sizeof(int), csize = sizeof(char), lsize = sizeof(long);

	printf("Size of a int is: %d\n", lsize);

	DataSize h = { 4096, 4096, n_Frames, 16 };

	clusterMap = (int *) malloc(h.nx * h.ny * isize);

	printf("There are %d frames. Starting to cluster on them..\n", h.nz);

	reduction_start = clock();
	for (z = 0; z < h.nz; z++)
	{
		int label = 0;
		int l0, l1, x0, y0, x1, y1;

		threshold_image (frameBuffer, darkBuffer, epsilon_s, h, clusterMap);
		
		/*
		//Label the clusters as -1, nonclusters 0
		for (y = 0; y < h.ny; y++) {
			for (x = 0; x < h.nx; x++) {
				l = y * h.nx + x;
				tmp_f = frameBuffer[l];
				tmp_d = darkBuffer[l];
				thresh_xy = tmp_d + epsilon;
				clusterMap[l] = (tmp_f > thresh_xy) ? -1 : 0;
			}
		}
		*/

		printf("Background subtraction done.\n");
		
		//forest_fire_conn_comp_labelling (frameBuffer, darkBuffer, epsilon_s, h, clusterMap, x_coor, y_coor);
		
		/*
		//set values of cluster map boundary to zero
		for (x = 0; x < h.nx; x++)
		{
			y0 = 0;
			y1 = Y_LEN - 1;
			l0 = y0 * h.nx + x;
			l1 = y1 * h.nx + x;
			clusterMap[l0] = 0;
			clusterMap[l1] = 0;
		}
		for (y = 0; y < h.ny; y++)
		{
			x0 = 0;
			x1 = h.nx - 1;
			l0 = y * h.nx + x0;
			l1 = y * h.nx + x1;
			clusterMap[l0] = 0;
			clusterMap[l1] = 0;
		}

		double stat_max = 0., stat_area = 0., stat_total = 0.;
		//Identify and label clusters
		for (y = 0; y < h.ny; y++)
		{
			for (x = 0; x < h.nx; x++)
			{
				int l1, l2, l_curr;
				l1 = y * h.nx + x;
				if (clusterMap[l1] == -1)
				{
					int n;
					int n_elem = 0, areaCount = 0, sum = 0;
					int x_curr, y_curr, x2, y2, maxVal;
					int nb[4][2] = { { -1,0 },{ 1,0 },{ 0,-1 },{ 0,1 } };
					int n_max = 4;
					bd_x[n_elem] = x;
					bd_y[n_elem] = y;
					maxVal = 0.;
					label++;
					do {
						//Consume the last pixel in search list
						x_curr = bd_x[n_elem];
						y_curr = bd_y[n_elem];
						l_curr = y_curr * h.nx + x_curr;
						clusterMap[l_curr] = label;
						areaCount++;
						dark_subtracted_pixel_val = frameBuffer[l_curr] - (darkBuffer[l_curr] + epsilon);
						dark_subtracted_pixel_val = (dark_subtracted_pixel_val > 0) ? dark_subtracted_pixel_val : 0;
						sum += dark_subtracted_pixel_val;
						if (dark_subtracted_pixel_val > maxVal)
						{
							maxVal = dark_subtracted_pixel_val;
							x_coor[label] = x_curr;
							y_coor[label] = y_curr;
						}
						n_elem--;

						//Add pixels connected to this cluster to future search list
						for (n = 0; n < n_max; n++)
						{
							x2 = x_curr + nb[n][0];
							y2 = y_curr + nb[n][1];
							l2 = y2 * h.nx + x2;
							if (clusterMap[l2] == -1)
							{
								n_elem++;
								bd_x[n_elem] = x2;
								bd_y[n_elem] = y2;
								clusterMap[l2] = label;
							}
						}

					} while (n_elem >= 0);

					max[label] = maxVal;
				}
			}
		}
		
		*/

		centroidImage = (unsigned char*)calloc((h.nx * h.ny) / 8, 1);
		for (x = 1; x <= label; x++) {
			l = y_coor[x] * h.nx + x_coor[x];
			SetBit(centroidImage, l);
		}
		printf("No. of Labels: %d\n", label);
		reduction_end = clock();
		reduction_time = (reduction_end - reduction_start) * 1000.0 / CLOCKS_PER_SEC;
		
		
		compression_start = clock();
		
		unsigned long nLenOrig = (h.nx * h.ny) / 8;
		unsigned long nCompressedSize;
		unsigned char* compressedCentroidImage;
		gzip_compress_stream (centroidImage, nLenOrig, mode, &nCompressedSize, &compressedCentroidImage);
		
		/*
		int nLenOrig = (h.nx * h.ny) / 8;
		printf("nLenOrig: %d\n", nLenOrig);
		int nLenDst = GetMaxCompressedLen(nLenOrig);

		//unsigned char* pbDst = new unsigned char[nLenDst];  // alloc dest buffer
		unsigned char* pbDst = (unsigned char*)calloc(nLenDst, 1)

		z_stream zInfo = { 0 };
		zInfo.total_in = zInfo.avail_in = nLenOrig;
		zInfo.total_out = zInfo.avail_out = nLenDst;
		zInfo.next_in = (Bytef *)centroidImage;
		zInfo.next_out = (Bytef *)pbDst;

		int nErr, nRet = -1;
		if (mode == 0) {
			nErr = deflateInit(&zInfo, Z_BEST_SPEED); // zlib function
		} else if(mode == 1) {
			nErr = deflateInit(&zInfo, Z_BEST_COMPRESSION); // zlib function
		}
		if (nErr == Z_OK) {
			nErr = deflate(&zInfo, Z_FINISH);              // zlib function
			if (nErr == Z_STREAM_END) {
				nRet = zInfo.total_out;
			}
		}
		deflateEnd(&zInfo);    // zlib function
		*/

		compression_end = clock();
		compression_time = (compression_end - compression_start) * 1000.0 / CLOCKS_PER_SEC;

		printf("Compression done.\n");
		//printf("Max Compressed Size: %d\n", nLenDst);
		printf("Compressed Size: %d\n", nCompressedSize);
		printf("Compression time for frame %d: %f ms\n", z, compression_time);

		compress_times_sizes[0] = reduction_time;
		compress_times_sizes[1] = compression_time;
		compress_times_sizes[2] = nLenOrig*1.0;
		compress_times_sizes[3] = nCompressedSize*1.0;
	}

	free(clusterMap);
	free(centroidImage);
	free(x_coor);
	free(y_coor);
	//free(max);
	

}

void do_L4_Reduce_Compress_Linear_Indices (unsigned short* frameBuffer, unsigned short* darkBuffer, int epsilon, int n_Frames, float* compress_times_sizes, int mode) {

	int error;
	int x, y, z;
	unsigned int l;

	unsigned short dark_subtracted_pixel_val;

	int *clusterMap;
	unsigned char *packed_indices;
	
	int tmp_f, tmp_d, thresh_xy;

	clock_t reduction_start, reduction_end, compression_start, compression_end;
	float reduction_time, compression_time;

	float *x_coor, *y_coor, *max;
	int *bd_x, *bd_y;

	x_coor = (float *) malloc(MAX_CLUSTER * sizeof(*x_coor));
	y_coor = (float *) malloc(MAX_CLUSTER * sizeof(*y_coor));
	max = (float *) malloc(MAX_CLUSTER * sizeof(*max));
	bd_x = (int *) malloc(MAX_ELEM * sizeof(*bd_x));
	bd_y = (int *) malloc(MAX_ELEM * sizeof(*bd_y));

	int fsize = sizeof(float), isize = sizeof(int), csize = sizeof(char), lsize = sizeof(long);

	printf("Size of a int is: %d\n", lsize);

	DataSize h = { 4096, 4096, n_Frames, 16 };

	clusterMap = (int *) malloc(h.nx * h.ny * isize);

	printf("There are %d frames. Starting to cluster on them..\n", h.nz);

	reduction_start = clock();
	for (z = 0; z < h.nz; z++)
	{
		int label = 0;
		int l0, l1, x0, y0, x1, y1;

		//Label the clusters as -1, nonclusters 0
		for (y = 0; y < h.ny; y++) {
			for (x = 0; x < h.nx; x++) {
				l = y * h.nx + x;
				tmp_f = frameBuffer[l];
				tmp_d = darkBuffer[l];
				thresh_xy = tmp_d + epsilon;
				clusterMap[l] = (tmp_f > thresh_xy) ? -1 : 0;
			}
		}

		printf("Background subtraction done.\n");

		//set values of cluster map boundary to zero
		for (x = 0; x < h.nx; x++)
		{
			y0 = 0;
			y1 = Y_LEN - 1;
			l0 = y0 * h.nx + x;
			l1 = y1 * h.nx + x;
			clusterMap[l0] = 0;
			clusterMap[l1] = 0;
		}
		for (y = 0; y < h.ny; y++)
		{
			x0 = 0;
			x1 = h.nx - 1;
			l0 = y * h.nx + x0;
			l1 = y * h.nx + x1;
			clusterMap[l0] = 0;
			clusterMap[l1] = 0;
		}

		double stat_max = 0., stat_area = 0., stat_total = 0.;
		//Identify and label clusters
		for (y = 0; y < h.ny; y++)
		{
			for (x = 0; x < h.nx; x++)
			{
				int l1, l2, l_curr;
				l1 = y * h.nx + x;
				if (clusterMap[l1] == -1)
				{
					int n;
					int n_elem = 0, areaCount = 0, sum = 0;
					int x_curr, y_curr, x2, y2, maxVal;
					int nb[4][2] = { { -1,0 },{ 1,0 },{ 0,-1 },{ 0,1 } };
					int n_max = 4;
					bd_x[n_elem] = x;
					bd_y[n_elem] = y;
					maxVal = 0.;
					label++;
					do {
						//Consume the last pixel in search list
						x_curr = bd_x[n_elem];
						y_curr = bd_y[n_elem];
						l_curr = y_curr * h.nx + x_curr;
						clusterMap[l_curr] = label;
						areaCount++;
						dark_subtracted_pixel_val = frameBuffer[l_curr] - (darkBuffer[l_curr] + epsilon);
						dark_subtracted_pixel_val = (dark_subtracted_pixel_val > 0) ? dark_subtracted_pixel_val : 0;
						sum += dark_subtracted_pixel_val;
						if (dark_subtracted_pixel_val > maxVal)
						{
							maxVal = dark_subtracted_pixel_val;
							x_coor[label] = x_curr;
							y_coor[label] = y_curr;
						}
						n_elem--;

						//Add pixels connected to this cluster to future search list
						for (n = 0; n < n_max; n++)
						{
							x2 = x_curr + nb[n][0];
							y2 = y_curr + nb[n][1];
							l2 = y2 * h.nx + x2;
							if (clusterMap[l2] == -1)
							{
								n_elem++;
								bd_x[n_elem] = x2;
								bd_y[n_elem] = y2;
								clusterMap[l2] = label;
							}
						}

					} while (n_elem >= 0);

					max[label] = maxVal;
				}
			}
		}

		packed_indices = (unsigned char *) calloc(label*3, 1);
		for (x = 1; x <= label; x++) {
			l = y_coor[x] * h.nx + x_coor[x];
			packed_indices[x*3] = l;
			packed_indices[x*3+1] = l/256;
			packed_indices[x*3+2] = l/65536;
		}
		printf("No. of Labels: %d\n", label);
		reduction_end = clock();
		reduction_time = (reduction_end - reduction_start) * 1000.0 / CLOCKS_PER_SEC;
		
		
		compression_start = clock();
		/*
		int nLenOrig = label*3;
		printf("nLenOrig: %d\n", nLenOrig);
		int nLenDst = GetMaxCompressedLen(nLenOrig);

		//unsigned char* pbDst = new unsigned char[nLenDst];  // alloc dest buffer
		unsigned char* pbDst = (unsigned char*)calloc(nLenDst, 1)

		z_stream zInfo = { 0 };
		zInfo.total_in = zInfo.avail_in = nLenOrig;
		zInfo.total_out = zInfo.avail_out = nLenDst;
		zInfo.next_in = (Bytef *)packed_indices;
		zInfo.next_out = (Bytef *)pbDst;

		int nErr, nRet = -1;
		if (mode == 0) {
			nErr = deflateInit(&zInfo, Z_BEST_SPEED); // zlib function
		} else if(mode == 1) {
			nErr = deflateInit(&zInfo, Z_BEST_COMPRESSION); // zlib function
		}
		if (nErr == Z_OK) {
			nErr = deflate(&zInfo, Z_FINISH);              // zlib function
			if (nErr == Z_STREAM_END) {
				nRet = zInfo.total_out;
			}
		}
		deflateEnd(&zInfo);    // zlib function
		*/

		compression_end = clock();
		compression_time = (compression_end - compression_start) * 1000.0 / CLOCKS_PER_SEC;

		printf("Compression done.\n");
		//printf("Max Compressed Size: %d\n", nLenDst);
		//printf("Compressed Size: %d\n", nRet);
		printf("Compression time for frame %d: %f ms\n", z, compression_time);

		compress_times_sizes[0] = reduction_time;
		compress_times_sizes[1] = compression_time;
		//compress_times_sizes[2] = nLenOrig*1.0;
		//compress_times_sizes[3] = nRet*1.0;
	}

	free(clusterMap);
	free(packed_indices);
	free(x_coor);
	free(y_coor);
	free(max);
	free(bd_x);
	free(bd_y);

}

void do_L4_Reduce_Compress_RLE (unsigned short* frameBuffer, unsigned short* darkBuffer, int epsilon, int n_Frames, float* compress_times_sizes, int mode) {

	int error;
	int x, y, z;
	unsigned int l;

	unsigned short dark_subtracted_pixel_val;

	int *clusterMap;
	unsigned char *centroidImage;
	unsigned char* rle;
	
	int tmp_f, tmp_d, thresh_xy;

	clock_t reduction_start, reduction_end, compression_start, compression_end;
	float reduction_time, compression_time;

	float *x_coor, *y_coor, *max;
	int *bd_x, *bd_y;

	x_coor = (float *) malloc(MAX_CLUSTER * sizeof(*x_coor));
	y_coor = (float *) malloc(MAX_CLUSTER * sizeof(*y_coor));
	max = (float *) malloc(MAX_CLUSTER * sizeof(*max));
	bd_x = (int *) malloc(MAX_ELEM * sizeof(*bd_x));
	bd_y = (int *) malloc(MAX_ELEM * sizeof(*bd_y));

	int fsize = sizeof(float), isize = sizeof(int), csize = sizeof(char), lsize = sizeof(long);

	printf("Size of a int is: %d\n", lsize);

	DataSize h = { 4096, 4096, n_Frames, 16 };

	clusterMap = (int *) malloc(h.nx * h.ny * isize);

	printf("There are %d frames. Starting to cluster on them..\n", h.nz);

	reduction_start = clock();
	for (z = 0; z < h.nz; z++)
	{
		int label = 0;
		int l0, l1, x0, y0, x1, y1;

		//Label the clusters as -1, nonclusters 0
		for (y = 0; y < h.ny; y++) {
			for (x = 0; x < h.nx; x++) {
				l = y * h.nx + x;
				tmp_f = frameBuffer[l];
				tmp_d = darkBuffer[l];
				thresh_xy = tmp_d + epsilon;
				clusterMap[l] = (tmp_f > thresh_xy) ? -1 : 0;
			}
		}

		printf("Background subtraction done.\n");

		//set values of cluster map boundary to zero
		for (x = 0; x < h.nx; x++)
		{
			y0 = 0;
			y1 = Y_LEN - 1;
			l0 = y0 * h.nx + x;
			l1 = y1 * h.nx + x;
			clusterMap[l0] = 0;
			clusterMap[l1] = 0;
		}
		for (y = 0; y < h.ny; y++)
		{
			x0 = 0;
			x1 = h.nx - 1;
			l0 = y * h.nx + x0;
			l1 = y * h.nx + x1;
			clusterMap[l0] = 0;
			clusterMap[l1] = 0;
		}

		double stat_max = 0., stat_area = 0., stat_total = 0.;
		//Identify and label clusters
		for (y = 0; y < h.ny; y++)
		{
			for (x = 0; x < h.nx; x++)
			{
				int l1, l2, l_curr;
				l1 = y * h.nx + x;
				if (clusterMap[l1] == -1)
				{
					int n;
					int n_elem = 0, areaCount = 0, sum = 0;
					int x_curr, y_curr, x2, y2, maxVal;
					int nb[4][2] = { { -1,0 },{ 1,0 },{ 0,-1 },{ 0,1 } };
					int n_max = 4;
					bd_x[n_elem] = x;
					bd_y[n_elem] = y;
					maxVal = 0.;
					label++;
					do {
						//Consume the last pixel in search list
						x_curr = bd_x[n_elem];
						y_curr = bd_y[n_elem];
						l_curr = y_curr * h.nx + x_curr;
						clusterMap[l_curr] = label;
						areaCount++;
						dark_subtracted_pixel_val = frameBuffer[l_curr] - (darkBuffer[l_curr] + epsilon);
						dark_subtracted_pixel_val = (dark_subtracted_pixel_val > 0) ? dark_subtracted_pixel_val : 0;
						sum += dark_subtracted_pixel_val;
						if (dark_subtracted_pixel_val > maxVal)
						{
							maxVal = dark_subtracted_pixel_val;
							x_coor[label] = x_curr;
							y_coor[label] = y_curr;
						}
						n_elem--;

						//Add pixels connected to this cluster to future search list
						for (n = 0; n < n_max; n++)
						{
							x2 = x_curr + nb[n][0];
							y2 = y_curr + nb[n][1];
							l2 = y2 * h.nx + x2;
							if (clusterMap[l2] == -1)
							{
								n_elem++;
								bd_x[n_elem] = x2;
								bd_y[n_elem] = y2;
								clusterMap[l2] = label;
							}
						}

					} while (n_elem >= 0);

					max[label] = maxVal;
				}
			}
		}

		centroidImage = (unsigned char*)calloc((h.nx * h.ny) / 8, 1);
		for (x = 1; x <= label; x++) {
			l = y_coor[x] * h.nx + x_coor[x];
			SetBit(centroidImage, l);
		}
		
		rle = (unsigned char *) calloc((label+1)*3, 1);
		int count = 0;
		int run_length = 0;
		for (x = 1; x < h.nx * h.ny; x++) {
			if (CheckBit(centroidImage, x)) {
				rle[count*3] = run_length;
				rle[count*3+1] = run_length/256;
				rle[count*3+2] = run_length/65536;
				count++;
				run_length = 0;
			} else {
				run_length++;
			}
		}
		printf("No. of Labels: %d\n", label);
		reduction_end = clock();
		reduction_time = (reduction_end - reduction_start) * 1000.0 / CLOCKS_PER_SEC;
		
		
		compression_start = clock();
		
		/*
		int nLenOrig = label*3;
		printf("nLenOrig: %d\n", nLenOrig);
		int nLenDst = GetMaxCompressedLen(nLenOrig);

		//unsigned char* pbDst = new unsigned char[nLenDst];  // alloc dest buffer
		unsigned char* pbDst = (unsigned char*)calloc(nLenDst, 1)

		z_stream zInfo = { 0 };
		zInfo.total_in = zInfo.avail_in = nLenOrig;
		zInfo.total_out = zInfo.avail_out = nLenDst;
		zInfo.next_in = (Bytef *)rle;
		zInfo.next_out = (Bytef *)pbDst;

		int nErr, nRet = -1;
		if (mode == 0) {
			nErr = deflateInit(&zInfo, Z_BEST_SPEED); // zlib function
		} else if(mode == 1) {
			nErr = deflateInit(&zInfo, Z_BEST_COMPRESSION); // zlib function
		}
		if (nErr == Z_OK) {
			nErr = deflate(&zInfo, Z_FINISH);              // zlib function
			if (nErr == Z_STREAM_END) {
				nRet = zInfo.total_out;
			}
		}
		deflateEnd(&zInfo);    // zlib function
		*/

		compression_end = clock();
		compression_time = (compression_end - compression_start) * 1000.0 / CLOCKS_PER_SEC;

		printf("Compression done.\n");
		//printf("Max Compressed Size: %d\n", nLenDst);
		//printf("Compressed Size: %d\n", nRet);
		printf("Compression time for frame %d: %f ms\n", z, compression_time);

		compress_times_sizes[0] = reduction_time;
		compress_times_sizes[1] = compression_time;
		//compress_times_sizes[2] = nLenOrig*1.0;
		//compress_times_sizes[3] = nRet*1.0;
	}

	free(clusterMap);
	free(centroidImage);
	free(rle);
	free(x_coor);
	free(y_coor);
	free(max);
	free(bd_x);
	free(bd_y);

}

void load_simulated_data(const char* inFolder, const char* darkFile, unsigned short* darkBuffer, unsigned short* frameBuffer) {
	
	int error;
	int z;
	int frame_start_index;
	
	FILE *fp_indata, *fp_darkdata;
	
	// load dark data
	fp_darkdata = fopen(darkFile, "rb");
	error = fread(darkBuffer, DTYPE_USHORT, FRAMEXDIM*FRAMEYDIM, fp_darkdata);
	printf("Dark data loaded.\n");
	
	// load frames
	for (z = 0; z < NUMFRAMES; z++) {
		
		char filename[WORDLEN];
		char temp[WORDLEN];
		sprintf(temp, "/Simulated_Frame_%d", z);
		strcpy(filename, inFolder);
		strcat(filename, temp);
		strcat(filename, ".bin");
		printf("%s\n", filename);
		
		fp_indata = fopen(filename, "rb");
		
		frame_start_index = z*(FRAMEXDIM*FRAMEYDIM);
		error = fread(frameBuffer + frame_start_index, DTYPE_USHORT, FRAMEXDIM*FRAMEYDIM, fp_indata);
	}
	
	fclose(fp_indata);
	fclose(fp_darkdata);
	
	printf("Image data loaded.\n");

}

void load_dataset (const char* filename, unsigned short* buffer, int nFrames, long int nOffset, DataSize d) {
	
	int error;

	FILE *fp;
	
	// load dark data
	fp = fopen(filename, "rb");
	fseek (fp, nOffset, SEEK_SET);
	error = fread(buffer, DTYPE_USHORT, d.ny*d.nx*nFrames, fp);
	fclose(fp);
	
	printf("Image data loaded.\n");

}

void showFrame (unsigned short* buffer, DataSize d, int framenum) {
	
	int r = 1;
	int c = 1;
	
	FILE * fp;

	fp = fopen ("C_Frame.txt", "w+");
	
	long int offset = d.ny*d.nx*framenum;
	for (r = 0; r < d.ny; r++) {
		for (c = 0; c < d.nx; c++) {
			fprintf(fp, "%d ", buffer[offset + r*d.nx + c]);
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
}

void logResult (const char* testname, float* compress_times_sizes) {
	
	FILE *logfile = fopen("results.txt", "a");
	float compression_ratio = compress_times_sizes[3]/compress_times_sizes[2];
	fprintf(logfile, "%s %2f %2f %2f %2f %2f\n", testname, compress_times_sizes[0], compress_times_sizes[1], compress_times_sizes[2], compress_times_sizes[3], compression_ratio);
	fclose(logfile);
}

void logResultHeader (const char* expname, const char *fieldnames[], int nFields) {

	/*
	FILE *logfile = fopen("results.txt", "a");
	fprintf(logfile, "\n");
	fprintf(logfile, "Experiment Name: Location Information Representation Compressibilities Test\n");
	fprintf(logfile, "Compression Representation Dose_Rate Reduction_Time Compression_Time Original_Size Compressed_Size Compression_Ratios\n");
	fclose(logfile);
	*/
	
	
	FILE *logfile1 = fopen("results.txt", "a");
	if (logfile1 == NULL)
	{
		printf("Error opening log file!\n");
	}
	fprintf(logfile1, "\n");
	
	fprintf(logfile1, "Experiment Name : %s\n", expname);
	int i;
	for (i=0; i<nFields; i++) {
		fprintf(logfile1, "%s ", fieldnames[i]);
	}
	
	fprintf(logfile1, "\n");
	fclose(logfile1);
	
}

int main(int argc, char *argv[]) {
	
	/*
	DataSize d = { 4096, 512, 11, 16 };
	unsigned long int b_size = (unsigned long int)d.nz * d.nx * d.ny * 2;
	unsigned short* buffer = (unsigned short*)malloc(b_size);
	load_dataset ("/home/abhik/code/DataCompression/beam_blanker_binary_data/12-23-00.232.bin", buffer, d.nz, 0, d);
	
	showFrame (buffer, d, 5);
	free(buffer);
	return 0;
	*/
	
	unsigned short* darkBuffer = (unsigned short*)malloc(FRAMEXDIM * FRAMEYDIM * 2);
	
	unsigned long int buffer_size = (unsigned long int)NUMFRAMES * FRAMEXDIM * FRAMEYDIM * 2;
	unsigned short* frameBuffer = (unsigned short*)malloc(buffer_size);
	
	const char *fieldnames[8] = {"Compression", "Representation", "Dose_Rate", "Reduction_Time", "Compression_Time", "Original_Size", "Compressed_Size", "Compression_Ratios"};
	logResultHeader ("Location_Information_Representation_Compressibility_Test", fieldnames, 8);
	
	const char *dose_rates[4] = {"0.8", "1.6", "3.2", "6.4"};
	
	int i;
	for (i=0; i<4; i++) {
		
		const char* inFolder = concat("/home/abhik/code/DataCompression/simulated_images/", dose_rates[i]);
		const char* darkFile = concat(concat("/home/abhik/code/DataCompression/simulated_images/", dose_rates[i]), "/Simulated_Dark_Frame.bin");
		load_simulated_data(inFolder, darkFile, darkBuffer, frameBuffer);
				
		int epsilon = 0;
		
		float compress_times_sizes_best_speed[4] = {0};
		float compress_times_sizes_best_size[4] = {0};
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 0, 12, compress_times_sizes_best_speed);		// best speed
		logResult (concat("GZIP-S L1_d12 ", dose_rates[i]), compress_times_sizes_best_speed);
		
		free(darkBuffer);
		free(frameBuffer);
		
		do_L1_Decompress_Expand ("Part_1.rc1");
		
		return 0;
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 1, 12, compress_times_sizes_best_size);		// best compression
		logResult (concat("GZIP-C L1_d12 ", dose_rates[i]), compress_times_sizes_best_size);
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 0, 10, compress_times_sizes_best_speed);		// best speed
		logResult (concat("GZIP-S L1_d10 ", dose_rates[i]), compress_times_sizes_best_speed);
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 1, 10, compress_times_sizes_best_size);		// best compression
		logResult (concat("GZIP-C L1_d10 ", dose_rates[i]), compress_times_sizes_best_size);
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 0, 8, compress_times_sizes_best_speed);		// best speed
		logResult (concat("GZIP-S L1_d8 ", dose_rates[i]), compress_times_sizes_best_speed);
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 1, 8, compress_times_sizes_best_size);		// best compression
		logResult (concat("GZIP-C L1_d8 ", dose_rates[i]), compress_times_sizes_best_size);
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 0, 4, compress_times_sizes_best_speed);		// best speed
		logResult (concat("GZIP-S L1_d4 ", dose_rates[i]), compress_times_sizes_best_speed);
		
		do_L1_Reduce_Compress (frameBuffer, darkBuffer, epsilon, 1, 1, 4, compress_times_sizes_best_size);		// best compression
		logResult (concat("GZIP-C L1_d4 ", dose_rates[i]), compress_times_sizes_best_size);
		
		/*
		do_L4_Reduce_Compress_Binary_Image (frameBuffer, darkBuffer, epsilon, 1, compress_times_sizes_best_speed, 0);		// best speed
		logResult (concat("GZIP-S Binary_Image ", dose_rates[i]), compress_times_sizes_best_speed);
		
		do_L4_Reduce_Compress_Binary_Image (frameBuffer, darkBuffer, epsilon, 1, compress_times_sizes_best_size, 1);		// best compression
		logResult (concat("GZIP-C Binary_Image ", dose_rates[i]), compress_times_sizes_best_size);
		
		do_L4_Reduce_Compress_Linear_Indices(frameBuffer, darkBuffer, epsilon, 1, compress_times_sizes_best_speed, 0);		// best speed
		logResult (concat("GZIP-S Linear_Indices ", dose_rates[i]), compress_times_sizes_best_speed);
		
		do_L4_Reduce_Compress_Linear_Indices (frameBuffer, darkBuffer, epsilon, 1, compress_times_sizes_best_size, 1);		// best compression
		logResult (concat("GZIP-C Linear_Indices ", dose_rates[i]), compress_times_sizes_best_size);
		
		do_L4_Reduce_Compress_RLE (frameBuffer, darkBuffer, epsilon, 1, compress_times_sizes_best_speed, 0);		// best speed
		logResult (concat("GZIP-S RLE ", dose_rates[i]), compress_times_sizes_best_speed);
		
		do_L4_Reduce_Compress_RLE (frameBuffer, darkBuffer, epsilon, 1, compress_times_sizes_best_size, 1);		// best compression
		logResult (concat("GZIP-C RLE ", dose_rates[i]), compress_times_sizes_best_size);
		*/
		
		free(darkBuffer);
		free(frameBuffer);
	}
}



