


void getDarkMax (uint16_t *darkBuffer, DataSize h, uint16_t *darkFrame, uint8_t nThreads) {
	
	uint64_t row, col, fno, pixel_index;
	uint64_t linear_index;
	uint64_t frame_size = h.nx*h.ny;
	
	printf("h.nz = %lu\n", h.nz);
	
	/*
	uint32_t count = 0;
	uint16_t dark_frames[h.nz][h.ny][h.nx];
	for (fno = 0; fno < h.nz; fno++) {
		for (row = 0; row < h.ny; row++) {
			for (col = 0; col < h.nx; col++) {
				dark_frames[fno][row][col] = 1;
				//printf("fno = %"PRIu32", row = %"PRIu32", col = %"PRIu32", p = %"PRIu16"\n", fno, row, col, darkMax_2[pixel_index]);
			}
		}
	}
	*/

	/*
	FILE *fp = fopen ("darks_pixel_00.txt", "w+");
	for (fno = 0; fno < h.nz; fno++) {
		row = 0;
		col = 0;
		pixel_index = row * h.nx + col;
		linear_index = h.nx*h.ny*fno + h.nx*row + col;
		fprintf(fp, "%"PRIu16"\n", darkBuffer[linear_index]);
	}
	fclose(fp);
	*/

	/*
	uint16_t *darkMax_2 = (uint16_t *)calloc(4096*512, sizeof(uint16_t));
	for (fno = 0; fno < h.nz; fno++) {
		for (row = 0; row < h.ny; row++) {
			for (col = 0; col < h.nx; col++) {
				pixel_index = row * h.nx + col;
				linear_index = h.nx*h.ny*fno + h.nx*row + col;
				if (darkBuffer[linear_index] > darkMax_2[pixel_index]) {
					darkMax_2[pixel_index] = darkBuffer[linear_index];
					if (row == 0 && col == 0) {
						printf("Dark Pixel (Max 2) %"PRIu64", %"PRIu64", %"PRIu64" = %"PRIu16"\n", fno, row, col, darkMax_2[pixel_index]);
					}
				}
			}
		}
	}
	*/
	
	/*
	for (row = 0; row < 25; row++) {
		for (col = 0; col < 25; col++) {
			pixel_index = row * h.nx + col;
			printf("Dark Pixel (Max 2) %lu, %lu = %d\n", row, col, darkMax_2[pixel_index]);
		}
	}
	*/
	
	//#pragma omp parallel for num_threads(nThreads) collapse(2)
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			pixel_index = row * h.nx + col;
			for (fno = 0; fno < h.nz; fno++) {
				linear_index = frame_size*fno + pixel_index;
				if (darkFrame[pixel_index] < darkBuffer[linear_index]) {
					darkFrame[pixel_index] = darkBuffer[linear_index];
				}
			}
		}
	}
	// implicit openmp barrier
	
	/*
	FILE *fp_max = fopen ("darks_max_00.txt", "w+");
	for (fno = 0; fno < h.nz; fno++) {
		row = 0;
		col = 0;
		pixel_index = row * h.nx + col;
		fprintf(fp_max, "%"PRIu16"\n", darkFrame[pixel_index]);
	}
	fclose(fp_max);
	*/
	
	/*
	FILE *fp_max = fopen ("darks_max_10000.txt", "w+");
	for (row = 0; row < 50; row++) {
		for (col = 0; col < 50; col++) {
			pixel_index = row * h.nx + col;
			fprintf(fp_max, "%"PRIu16"\n", darkFrame[pixel_index]);
		}
	}
	fclose(fp_max);
	*/
	
	/*
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			pixel_index = row * h.nx + col;
			for (fno = 0; fno < h.nz; fno++) {
				linear_index = frame_size*fno + pixel_index;
				printf("%d, %d, %d: %hu\n", row, col, fno, darkBuffer[linear_index]);
			}
		}
	}
	*/
	
	/*
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			pixel_index = row * h.nx + col;
			printf("%d, %d: %hu\n", row, col, darkFrame[pixel_index]);
		}
	}
	*/
	
}

void loadBinaryData (const char* filename, uint16_t *buffer, uint32_t frameStart, uint32_t n_Frames, DataSize d) {
	
	// load data
	FILE *fp = fopen(filename, "rb");
	size_t seek_count = fseek(fp, d.nx*d.ny*(uint64_t)frameStart, SEEK_SET);
	
	uint64_t n_pixels_in_frame = d.nx*d.ny;
	uint64_t chunk_size_frames = 200;
	uint64_t chunk_size_pixels = chunk_size_frames * n_pixels_in_frame;
	uint64_t iter_c = floor((uint64_t)n_Frames / chunk_size_frames);
	uint64_t rem_sz	= (uint64_t)n_Frames % chunk_size_frames;
	uint64_t i;
	uint64_t read_count = 0;
	uint64_t n_elem = 0;
	
	recode_print("nx: %lu\n", d.nx);
	recode_print("ny: %lu\n", d.ny);
	recode_print("nf: %"PRIu64"\n", n_Frames);
	recode_print("iter_c: %"PRIu64"\n", iter_c);
	
	for (i = 0; i < iter_c; i++) {
		n_elem = i*chunk_size_pixels;
		read_count += fread(buffer + n_elem, sizeof(uint16_t), chunk_size_pixels, fp);
	}
	n_elem = i*chunk_size_pixels;
	read_count += fread(buffer + i*chunk_size_pixels, sizeof(uint16_t), n_pixels_in_frame*rem_sz, fp);
	fclose(fp);
	
	recode_print("Seek Count: %"PRIu64"\n", seek_count);
	recode_print("Read Count: %"PRIu64"\n", read_count);
	
	uint64_t expected_read_count = d.nx*d.ny*n_Frames;
	recode_print("Expected Read Count: %"PRIu64"\n", expected_read_count);
	
	if (read_count != expected_read_count) {
		printf("An unexpected error occurred.\nActual and expected read counts do not match. Exiting.");
		exit(0);
	}
}

void loadMRCData(const char* filename, uint16_t *buffer, uint32_t frameStart, uint32_t n_Frames, DataSize d) {
	MRCHeader *mrc_header = (MRCHeader *)malloc(sizeof(MRCHeader));
	//parseMRCHeader(filename, &mrc_header);
	// load data
	FILE *fp = fopen(filename, "rb");
	getMRCFrames(fp, buffer, frameStart, n_Frames, d);
}

uint64_t loadSEQData(const char* filename, uint16_t *buffer, uint32_t frameStart, uint32_t n_Frames) {
	FILE *fp = fopen(filename, "rb");
	if (!fp) {
		printf("Failed to open file: '%s' for reading.", filename);
		exit(1);
	}
	SEQHeader *seq_header = (SEQHeader *)malloc(sizeof(SEQHeader));
	getSEQHeader(fp, &seq_header);
	// load data
	uint64_t retval = getSEQFrames(0, fp, buffer, frameStart, n_Frames, seq_header);
	fclose(fp);
	return retval;
}

uint64_t loadData(uint8_t file_type, const char* filename, uint16_t *buffer, uint32_t frameStart, uint32_t n_Frames, DataSize d) {

	if (file_type == SOURCE_FILE_TYPE_BINARY) {
		loadBinaryData(filename, buffer, frameStart, n_Frames, d);
		return 0;
	}
	else if (file_type == SOURCE_FILE_TYPE_MRC) {
		loadMRCData(filename, buffer, frameStart, n_Frames, d);
		return 0;
	}
	else if (file_type == SOURCE_FILE_TYPE_SEQUENCE) {
		return loadSEQData(filename, buffer, frameStart, n_Frames);
	}

}

// Used in L4 and L2
float get_foreground_image (uint16_t *frameBuffer, 
							uint16_t *darkBuffer, 
							uint16_t epsilon, 
							DataSize h, 
							uint32_t *n_fg_pixels, 
							uint16_t *foregroundImage, 
							int8_t   *foregroundTernaryMap) {
	
	clock_t p_start = clock();
	
	uint16_t row, col;
	uint32_t linear_index;
	uint16_t tmp_f, tmp_d, thresh_xy, dark_subtracted_pixval;
	
	*n_fg_pixels = 0;
	
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			tmp_d = darkBuffer[linear_index];
			thresh_xy = tmp_d + epsilon;

			if (tmp_f > thresh_xy) {
				foregroundImage [linear_index] = tmp_f - thresh_xy;
				foregroundTernaryMap [linear_index] = -1;
				(*n_fg_pixels)++;
			} else {
				foregroundImage [linear_index] = 0;
				foregroundTernaryMap [linear_index] = 0;
			}
			
		}
	}
	
	//printf("No. of foreground pixels = %lu\n", *nFgPixels);

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}




