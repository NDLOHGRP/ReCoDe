

/*
this function is called by pyrecode_c.py
*/
int64_t decompressExpand_L1_Reduced_Compressed_Frame_Sparse(
	FILE* rc_fp, uint16_t nx, uint16_t ny, uint8_t bit_depth, 
	uint32_t n_compressed_bytes_in_binary_image, uint32_t n_compressed_bytes_in_pixvals, uint32_t n_bytes_in_packed_pixvals, uint32_t n_bytes_in_binary_image,
	uint8_t *compressedBinaryImage, uint8_t *deCompressedBinaryImage, uint8_t *compressedPixvals, uint8_t *deCompressedPixvals,
	uint64_t *pow2_lookup_table,
	uint16_t *frameBuffer,
	uint8_t is_intermediate_file,
	uint8_t get_frame_id
) {

	/*
	rc_fp is assumed to point to beginning of a frame. Appropriate seek is handled by callers: pyrecode_c.py
	*/
	size_t t;
	uint32_t frame_id;
	if (is_intermediate_file == 1) {
		// skip frame_id, nCompressedSize_BinaryImage, nCompressedSize_Pixvals, bytesRequiredForPacking
		if (get_frame_id) {
			t = fread(&frame_id, sizeof(uint32_t), 1, rc_fp);
			if (t != 1) {
				return -1;
			} else {
				return (int64_t)frame_id;
			}
		}
		t = fread(&n_compressed_bytes_in_binary_image, sizeof(uint32_t), 1, rc_fp);
		if (t != 1) {
			return -1;
		}
		t = fread(&n_compressed_bytes_in_pixvals, sizeof(uint32_t), 1, rc_fp);
		if (t != 1) {
			return -1;
		}
		t = fread(&n_bytes_in_packed_pixvals, sizeof(uint32_t), 1, rc_fp);
		if (t != 1) {
			return -1;
		}
	}

	/*
	uint8_t byte_depth = ceil((bit_depth*1.0) / 8.0);
	uint32_t n_pixels_in_frame = (uint32_t)nx * (uint32_t)ny;					// number of pixels in the image
	uint32_t n_bytes_in_image = ceil(n_pixels_in_frame * byte_depth);			// number of bytes needed to hold pixvals
	compressedBinaryImage = (uint8_t*)calloc(100000, sizeof(uint8_t));
	deCompressedBinaryImage = (uint8_t*)calloc(100000, sizeof(uint8_t));
	compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	deCompressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	*/
	
	/*
	printf("L1.h: Processing Frame %d \n", frame_id);
	printf("nx = %d\n", nx);
	printf("ny = %d\n", ny);
	printf("n_compressed_bytes_in_binary_image = %" PRIu32 "\n", n_compressed_bytes_in_binary_image);
	printf("n_compressed_bytes_in_pixvals = %" PRIu32 "\n", n_compressed_bytes_in_pixvals);
	printf("n_bytes_in_packed_pixvals = %" PRIu32 "\n", n_bytes_in_packed_pixvals);
	printf("current position (L1.h): %d\n", ftell(rc_fp));
	*/

	size_t result_1 = fread(compressedBinaryImage, sizeof(uint8_t), n_compressed_bytes_in_binary_image, rc_fp);
	if (result_1 != n_compressed_bytes_in_binary_image) {
		return -1;
	}
	
	/*
	printf("L1.h: read %zu\n", result_1);
	printf("current position (L1.h): %d\n", ftell(rc_fp));
	*/

	size_t result_2 = fread(compressedPixvals, sizeof(uint8_t), n_compressed_bytes_in_pixvals, rc_fp);
	if (result_2 != n_compressed_bytes_in_pixvals) {
		return -1;
	}
	
	/*
	printf("L1.h: read %zu\n", result_2);
	printf("current position (L1.h): %d\n", ftell(rc_fp));
	for (int i = 0; i < 10; i++) {
		printf("compressedBinaryImage[%d] = %" PRIu8 "\n", i, compressedBinaryImage[i]);
	}
	*/
	
	decompress_stream(-1, -1, compressedBinaryImage, deCompressedBinaryImage, n_compressed_bytes_in_binary_image, n_bytes_in_binary_image);
	decompress_stream(-1, -1, compressedPixvals, deCompressedPixvals, n_compressed_bytes_in_pixvals, n_bytes_in_packed_pixvals);

	uint16_t row, col;
	uint32_t linear_pixel_index, pixel_bit_index_frame;
	uint16_t extracted_pixval;
	uint8_t n;
	uint64_t n_fg_pixels = 0;

	/*
	for (int i = 0; i < (4096 * 512) / 8; i++) {
		if (deCompressedBinaryImage[i] > 0) {
			printf("deCompressedBinaryImage[%d] = %" PRIu8 "\n", i, deCompressedBinaryImage[i]);
		}
	}
	*/

	for (row = 0; row < ny; row++) {
		for (col = 0; col < nx; col++) {
			linear_pixel_index = row * nx + col;
			if (CheckBit(deCompressedBinaryImage, linear_pixel_index) > 0) {
				// unpack pixval
				extracted_pixval = 0;
				// Assumes LITTLE-ENDIAN Byte Order of pixval
				for (n = 0; n < bit_depth; n++) {
					pixel_bit_index_frame = n_fg_pixels*bit_depth + n;
					if (CheckBit(deCompressedPixvals, pixel_bit_index_frame) > 0) {
						extracted_pixval += pow2_lookup_table[n];
					}
				}
				frameBuffer[n_fg_pixels * 3] = row;
				frameBuffer[n_fg_pixels * 3 + 1] = col;
				frameBuffer[n_fg_pixels * 3 + 2] = extracted_pixval;

				/*==============DEBUG ONLY===============
				if (row == 2830 && col == 3173) {
					//printf("Row = %" PRIu16 ", Col = %" PRIu16 ", Value = %" PRIu16 ", foreground pixel = %" PRIu64 "\n", row, col, extracted_pixval, n_fg_pixels);
					printf("Row = %d, Col = %d, Value = %d, foreground pixel = %d\n", row, col, extracted_pixval, n_fg_pixels);
					for (n = 0; n < bit_depth; n++) {
						pixel_bit_index_frame = n_fg_pixels*bit_depth + n;
						if (CheckBit(deCompressedPixvals, pixel_bit_index_frame) > 0) {
							printf("1 ");
						} else {
							printf("0 ");
						}
						printf("\n");
					}
				}
				==============DEBUG ONLY===============*/

				n_fg_pixels++;
			}
		}
	}
	//printf("Decoded Frame with %d foreground pixels\n", n_fg_pixels);
	return (int64_t)n_fg_pixels;
}

/*
These function are provided for testing purposes only. No external pythonic calls to these function are available. 
Decompression to only sparse format is supported for external calls.
void decompressExpand_L1_Reduced_Only(FILE* fp, const char *out_fname, RCHeader *header) {}
void decompressExpand_L1_Reduced_Only_Sparse(FILE* fp, const char *out_fname, RCHeader *header) {}
*/


/*
Called by the python function _bit_pack_pixel_intensities
a similar function with scaling is declared above: scale_and_pack_pixvals
*/
float bit_pack_pixel_intensities (
						uint64_t sz_packedPixval, 
						uint32_t n_fg_pixels,
						uint8_t  bit_depth, 
						uint16_t *pixvals, 
						uint8_t  *packedPixvals) {
	
	clock_t p_start = clock();

	int p, n, linear_index, nth_bit_of_pth_pixval;
	
	// setting packedPixvals to 0 in a for loop is faster than using ClearBit
	for (p = 0; p < sz_packedPixval; p++) {
		packedPixvals[p] = 0;
	}
	
	for (p=0; p<n_fg_pixels; p++) {
		
		// Assumes LITTLE-ENDIAN Byte Order
		for (n=0; n<bit_depth; n++) {
			nth_bit_of_pth_pixval = (pixvals[p] & ( 1 << n )) >> n;
			if (nth_bit_of_pth_pixval != 0) {
				linear_index = p*bit_depth + n;
				SetBit(packedPixvals, linear_index);
			}
		}
		
		//printf("%f, %hu, %hu, %hu, %hu\n", pixval_01, scaled_pixval, pixvals[p], data_min, data_max);
	}
	
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
	
}
