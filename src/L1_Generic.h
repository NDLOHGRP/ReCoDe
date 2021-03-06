


// Used in L1

#define BINARIZE_AND_GET_PIXVALS_TEMPLATE( TYPE )          							\
  float binarize_and_get_pixvals_##TYPE( TYPE *frameBuffer, 						\
  										 TYPE *darkBuffer, 							\
  										 TYPE epsilon, 								\
  										 DataSize h, 								\
  										 uint8_t *binaryImage, 						\
  										 TYPE *pixvals, 							\
  										 uint32_t *n_fg_pixels, 					\
  										 TYPE *data_min, 							\
  										 TYPE *data_max, 							\
  										 uint8_t process_id) 						\
  {                                  												\
  	clock_t p_start = clock();														\
    uint16_t row, col;																\
	uint32_t linear_index;															\
	TYPE	 tmp_f, tmp_d, thresh_xy, dark_subtracted_pixval;						\
																					\
	*n_fg_pixels = 0;																\
	*data_min = 65535;																\
	*data_max = 0;																	\
																					\
	for (row = 0; row < h.ny; row++) {												\
		for (col = 0; col < h.nx; col++) {											\
																					\
			linear_index = row * h.nx + col;										\
			tmp_f = frameBuffer[linear_index];										\
			tmp_d = darkBuffer[linear_index];										\
			thresh_xy = tmp_d + epsilon;											\
																					\
			if (tmp_f > thresh_xy) {												\
				SetBit(binaryImage, linear_index);									\
				dark_subtracted_pixval = tmp_f - thresh_xy;							\
				pixvals[(*n_fg_pixels)++] = dark_subtracted_pixval;					\
				if (dark_subtracted_pixval > *data_max) {							\
					*data_max = dark_subtracted_pixval;								\
				}																	\
				if (dark_subtracted_pixval < *data_min) {							\
					*data_min = dark_subtracted_pixval;								\
				}																	\
			} else {																\
				ClearBit(binaryImage, linear_index);								\
				*data_min = 0;														\
			}																		\
		}																			\
	}																				\
																					\
	printf("No. of foreground pixels = %lu\n", *n_fg_pixels);						\
																					\
	clock_t p_end = clock();														\
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;				\
																					\
    return process_time;          													\
  }


BINARIZE_AND_GET_PIXVALS_TEMPLATE (uint8_t)
BINARIZE_AND_GET_PIXVALS_TEMPLATE (uint16_t)

/*
binarize_and_get_pixvals_uint8_t (uint8_t *frameBuffer, 
								  uint8_t *darkBuffer, 
								  uint8_t epsilon, 
								  DataSize h, 
								  uint8_t  *binaryImage, 
								  uint8_t *pixvals, 
								  uint32_t *n_fg_pixels, 
								  uint8_t *data_min, 
								  uint8_t *data_max,
								  uint8_t	 process_id) {
	
	clock_t p_start = clock();
	
	uint16_t row, col;
	uint32_t linear_index;
	uint16_t tmp_f, tmp_d, thresh_xy, dark_subtracted_pixval;
	
	*n_fg_pixels = 0;
	*data_min = 65535;
	*data_max = 0;
	//printf("RCT %d: 1.2.1.0\n", process_id);
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			//printf("RCT %d: 1.2.1.1\n", process_id);
			tmp_d = darkBuffer[linear_index];
			//printf("RCT %d: 1.2.1.2.\n", process_id);
			thresh_xy = tmp_d + epsilon;

			//printf("%hu, %hu, %hu, %hu, %hu\n", row, col, tmp_f, tmp_d, thresh_xy);
			
			if (tmp_f > thresh_xy) {
				SetBit(binaryImage, linear_index);
				dark_subtracted_pixval = tmp_f - thresh_xy;
				pixvals[(*n_fg_pixels)++] = dark_subtracted_pixval;
				if (dark_subtracted_pixval > *data_max) {
					*data_max = dark_subtracted_pixval;
					//printf("MAX: %hu\n", *data_max);
				}
				if (dark_subtracted_pixval < *data_min) {
					*data_min = dark_subtracted_pixval;
				}
				//printf("row: %hu, col: %hu, val = %hu, dark = %hu, count = %hu, index = %"PRIu32", val = %hu\n", 
				//		row, col, tmp_f, thresh_xy, (*n_fg_pixels), linear_index, dark_subtracted_pixval);
			} else {
				ClearBit(binaryImage, linear_index);
				*data_min = 0;
				//printf("row: %hu, col: %hu, val = %hu, dark = %hu, count = %hu\n", 
				//		row, col, tmp_f, thresh_xy, (*n_fg_pixels));
			}
			//printf("RCT %d: 1.2.1.3\n", process_id);
		}
	}
	
	printf("No. of foreground pixels = %lu\n", *n_fg_pixels);

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	//return process_time;
}


binarize_and_get_pixvals_uint16_t (uint8_t	 process_id,
								uint16_t *frameBuffer, 
								uint16_t *darkBuffer, 
								uint16_t epsilon, 
								DataSize h, 
								uint8_t  *binaryImage, 
								uint16_t *pixvals, 
								uint32_t *n_fg_pixels, 
								uint16_t *data_min, 
								uint16_t *data_max) {
	
	clock_t p_start = clock();
	
	uint16_t row, col;
	uint32_t linear_index;
	uint16_t tmp_f, tmp_d, thresh_xy, dark_subtracted_pixval;
	
	*n_fg_pixels = 0;
	*data_min = 65535;
	*data_max = 0;
	//printf("RCT %d: 1.2.1.0\n", process_id);
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			
			linear_index = row * h.nx + col;
			tmp_f = frameBuffer[linear_index];
			//printf("RCT %d: 1.2.1.1\n", process_id);
			tmp_d = darkBuffer[linear_index];
			//printf("RCT %d: 1.2.1.2.\n", process_id);
			thresh_xy = tmp_d + epsilon;

			//printf("%hu, %hu, %hu, %hu, %hu\n", row, col, tmp_f, tmp_d, thresh_xy);
			
			if (tmp_f > thresh_xy) {
				SetBit(binaryImage, linear_index);
				dark_subtracted_pixval = tmp_f - thresh_xy;
				pixvals[(*n_fg_pixels)++] = dark_subtracted_pixval;
				if (dark_subtracted_pixval > *data_max) {
					*data_max = dark_subtracted_pixval;
					//printf("MAX: %hu\n", *data_max);
				}
				if (dark_subtracted_pixval < *data_min) {
					*data_min = dark_subtracted_pixval;
				}
				//printf("row: %hu, col: %hu, val = %hu, dark = %hu, count = %hu, index = %"PRIu32", val = %hu\n", 
				//		row, col, tmp_f, thresh_xy, (*n_fg_pixels), linear_index, dark_subtracted_pixval);
			} else {
				ClearBit(binaryImage, linear_index);
				*data_min = 0;
				//printf("row: %hu, col: %hu, val = %hu, dark = %hu, count = %hu\n", 
				//		row, col, tmp_f, thresh_xy, (*n_fg_pixels));
			}
			//printf("RCT %d: 1.2.1.3\n", process_id);
		}
	}
	
	printf("No. of foreground pixels = %lu\n", *n_fg_pixels);

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	//return process_time;
}


#define binarize_and_get_pixvals (frameBuffer, darkBuffer, epsilon, h, binaryImage, pixvals, n_fg_pixels, data_min, data_max, process_id)		\
	_Generic((frameBuffer),																														\
		uint8_t:  binarize_and_get_pixvals_8,			\
		uint16_t: binarize_and_get_pixvals_16,			\
		uint32_t: binarize_and_get_pixvals_32,			\
		uint64_t: binarize_and_get_pixvals_64			\
			)(frameBuffer, darkBuffer, epsilon, h, binaryImage, pixvals, n_fg_pixels, data_min, data_max, process_id)
*/

float scale_and_pack_pixvals (	uint8_t	 process_id,
								uint16_t *pixvals, 
								uint32_t n_fg_pixels, 
								uint16_t data_min, 
								uint16_t data_max, 
								uint8_t  bit_depth, 
								uint8_t  *packedPixvals) {
	
	clock_t p_start = clock();

	int p, n, linear_index, nth_bit_of_pth_pixval;
	unsigned short scaled_pixval;
	double pixval_01;
	
	unsigned short MAX_VAL = pow(2, bit_depth) - 1;
	
	int doScaling = 0;
	if (data_max > MAX_VAL) {
		doScaling = 1;
	}
	
	for (p=0; p<n_fg_pixels; p++) {
		
		if (doScaling) {
			pixval_01 = ((pixvals[p] - data_min)*1.0)/((data_max - data_min)*1.0);
			scaled_pixval = (unsigned short)ceil(pixval_01*MAX_VAL);
		} else {
			scaled_pixval = pixvals[p];
		}
		
		// Assumes LITTLE-ENDIAN Byte Order
		for (n=0; n<bit_depth; n++) {
			nth_bit_of_pth_pixval = (scaled_pixval & ( 1 << n )) >> n;
			if (nth_bit_of_pth_pixval == 1) {
				linear_index = p*bit_depth + n;
				SetBit(packedPixvals, linear_index);
			}
		}
		
		//printf("%f, %hu, %hu, %hu, %hu\n", pixval_01, scaled_pixval, pixvals[p], data_min, data_max);
	}
	
	//printf("MAX_VAL for given bit_depth: %hu\n", MAX_VAL);
	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
	
}



float reduceFrame_L1 (  uint8_t	 process_id,
						uint16_t *frameBuffer, 
						uint16_t *darkFrame, 
						uint32_t z, 
						DataSize h, 
						uint8_t  epsilon_s,
						uint8_t  bit_depth,
						uint16_t *pixvals, 
						uint8_t  *binaryImage,
						uint8_t  *packedPixvals,
						uint32_t *n_bytes_in_packed_pixvals) {
	
	
	/*
	========================================================================================== 
	Reduce Frame to Level 1
	========================================================================================== 
	*/
	
	uint32_t n_fg_pixels 		= 0;
	uint16_t data_min 	 		= 65535;
	uint16_t data_max 	 		= 0;
	uint32_t frame_start_index	= z * h.nx * h.ny;
		
	//printf("RCT %d: frame_start_index = %" PRIu32 "\n", process_id, frame_start_index);
	float thresh_time 			= binarize_and_get_pixvals_uint16_t (frameBuffer + frame_start_index, darkFrame, epsilon_s, h, 
															binaryImage, pixvals, &n_fg_pixels, &data_min, &data_max, process_id);
	//printf("RCT %d: 1.2.1\n", process_id);
	*n_bytes_in_packed_pixvals 	= ceil((n_fg_pixels * (bit_depth))/8);
	float packing_time 			= scale_and_pack_pixvals (process_id, pixvals, n_fg_pixels, data_min, data_max, bit_depth, packedPixvals);
	//printf("RCT %d: 1.2.2.\n", process_id);
	float reduction_time 		= thresh_time + packing_time;
	
	return reduction_time;
}



float compressFrame_L1 (uint8_t	 process_id,
						uint8_t  *binaryImage,
						uint8_t  *packedPixvals,
						uint32_t *n_bytes_in_packed_pixvals,
						DataSize h, 
						uint8_t  compression_scheme,
						uint8_t  compression_level,
						uint8_t  *compressedBinaryImage,
						uint8_t  *compressedPackedPixvals,
						uint32_t *n_compressed_bytes_1,
						uint32_t *n_compressed_bytes_2 ) {
	
	/*
	========================================================================================== 
	Compress
	========================================================================================== 
	*/
	uint32_t n_bytes_in_binary_image = ceil((h.nx * h.ny) / 8.0);
	
	float 	 compression_time_1 = compress_stream (compression_scheme, compression_level, binaryImage,   n_bytes_in_binary_image, n_compressed_bytes_1, compressedBinaryImage);
	float 	 compression_time_2 = compress_stream (compression_scheme, compression_level, packedPixvals, *n_bytes_in_packed_pixvals, n_compressed_bytes_2, compressedPackedPixvals);
	
	float 	 compression_time 	= compression_time_1 + compression_time_2;
	
	return compression_time;
}


float decompressFrame_L1 (	uint8_t	 process_id,
							DataSize h, 
							uint8_t  compression_scheme,
							uint8_t  compression_level,
							uint8_t  *compressedBinaryImage,
							uint8_t  *compressedPackedPixvals,
							uint32_t *n_compressed_bytes_1,
							uint32_t *n_compressed_bytes_2 ) {
	
	/*
	========================================================================================== 
	Decompress
	========================================================================================== 
	*/
	uint32_t n_pixels_in_frame			= h.nx * h.ny;														// number of pixels in the image
	uint32_t n_bytes_in_binary_image 	= ceil (n_pixels_in_frame / 8.0);									// number of bytes needed to pack binary image

	uint8_t  *temp 				 		= (uint8_t*)malloc (n_pixels_in_frame * sizeof(uint8_t) * 1.5);

	float 	 decompression_time_1 		= decompress_stream (compression_scheme, compression_level, compressedBinaryImage, temp, *n_compressed_bytes_1, -1);
	float 	 decompression_time_2 		= decompress_stream (compression_scheme, compression_level, compressedPackedPixvals, temp, *n_compressed_bytes_1, -1);
	float 	 decompression_time 		= decompression_time_1 + decompression_time_2;
	
	return decompression_time;
}


void reduceCompressFrame_L1 (	uint8_t	 process_id,
								uint16_t *frameBuffer, 
								uint16_t *darkFrame, 
								uint32_t frame_id,					// absolute frame index
								uint32_t z, 						// frame index relative to start of this thread
								DataSize h, 
								uint8_t  rc_operation_mode,
								uint8_t  epsilon_s,
								uint8_t  bit_depth,
								uint8_t  compression_scheme,
								uint8_t  compression_level,
								uint16_t *pixvals, 
								uint8_t  *binaryImage, 
								uint8_t  *packedPixvals, 
								uint8_t  *compressedBinaryImage,
								uint8_t  *compressedPackedPixvals,
								uint32_t *n_bytes_in_packed_pixvals, 
								uint32_t *n_compressed_bytes_1, 
								uint32_t *n_compressed_bytes_2,
								float 	 *run_metrics,
								uint8_t	 copy_compressed_frame_to_return_buffer,
								uint8_t  *compressed_frame,
								uint32_t *compressed_frame_length) {

		//printf("RCT %d: compressed_frame[0] = %d\n", process_id, compressed_frame[0]);
								
		//printf("RCT %d: 1.1.\n", process_id);
		float reduction_time  	= reduceFrame_L1 (process_id, frameBuffer, darkFrame, z, h, 
												  epsilon_s, bit_depth, 
												  pixvals, binaryImage, packedPixvals, n_bytes_in_packed_pixvals);
		//printf("RCT %d: 1.2.\n", process_id);
		float compression_time  	= 0.0;
		float decompression_time  	= 0.0;

		if (rc_operation_mode == RC_MODE_REDUCE_COMPRESS) {
												  
			compression_time  	= compressFrame_L1 (process_id, binaryImage, packedPixvals, n_bytes_in_packed_pixvals, h,
													compression_scheme, compression_level,
						  							compressedBinaryImage, compressedPackedPixvals, 
						  							n_compressed_bytes_1, n_compressed_bytes_2);
			//printf("RCT %d: 1.3.\n", process_id);

			
			decompression_time	= decompressFrame_L1 (	process_id, h,
														compression_scheme, compression_level,
						  								compressedBinaryImage, compressedPackedPixvals, 
						  								n_compressed_bytes_1, n_compressed_bytes_2);
			
		}

		printf("RCT %"PRIu8": Frame ID = %"PRIu32", n_bytes_in_packed_pixvals = %"PRIu32"\n", process_id, frame_id, *n_bytes_in_packed_pixvals);

		uint32_t n_bytes_in_binary_image = ceil((h.nx * h.ny) / 8.0);

		run_metrics[0] += reduction_time;
		run_metrics[1] += compression_time;
		run_metrics[2] += (n_bytes_in_binary_image + *n_bytes_in_packed_pixvals);
		run_metrics[3] += (*n_compressed_bytes_1 + *n_compressed_bytes_2);
		run_metrics[4] += decompression_time;
		
		
		/*
		========================================================================================== 
		Copy compressed data into return buffer if requested. Used in on-the-fly compression.
		========================================================================================== 
		*/
		if (copy_compressed_frame_to_return_buffer == 1) {
			
			if (rc_operation_mode 			== RC_MODE_REDUCE_ONLY) {
			
				*compressed_frame_length 	= sizeof(uint8_t) * (n_bytes_in_binary_image) + sizeof(uint8_t) * (*n_bytes_in_packed_pixvals) + 2 * sizeof(uint32_t);

				memcpy(compressed_frame, &frame_id, 					sizeof(uint32_t));
				memcpy(compressed_frame, n_bytes_in_packed_pixvals, 	sizeof(uint32_t));
				memcpy(compressed_frame, binaryImage, 					sizeof(uint8_t) * n_bytes_in_binary_image);
				memcpy(compressed_frame, packedPixvals, 				sizeof(uint8_t) * (*n_bytes_in_packed_pixvals));
				
				///printf("Copying to return buffer in Reduce Only Mode.\n");
			
			} else if (rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {
				
				*compressed_frame_length 	= sizeof(uint8_t) * (*n_compressed_bytes_1) + sizeof(uint8_t) * (*n_compressed_bytes_2) + 4 * sizeof(uint32_t);
				
				memcpy(compressed_frame, &frame_id, 					sizeof(uint32_t));
				memcpy(compressed_frame, n_compressed_bytes_1, 			sizeof(uint32_t));
				memcpy(compressed_frame, n_compressed_bytes_2, 			sizeof(uint32_t));
				memcpy(compressed_frame, n_bytes_in_packed_pixvals, 	sizeof(uint32_t));
				memcpy(compressed_frame, compressedBinaryImage, 		sizeof(uint8_t)*(*n_compressed_bytes_1));
				memcpy(compressed_frame, compressedPackedPixvals, 		sizeof(uint8_t)*(*n_compressed_bytes_2));
				
				///printf("Copying to return buffer in Reduce Compress Mode.\n");
			}
		}
		
		return;
	
}



char* reduceCompress_L1 (uint8_t	 process_id, 
						 const char  *original_filename, 
						 uint16_t    *frameBuffer, 				// points to the start of this thread's chunk
						 uint16_t    *darkFrame, 
						 uint32_t    frame_start_index, 
						 DataSize    h, 
						 InputParams *input_params,
						 float 		 *compression_time ) {
	
	
		
		uint32_t n_pixels_in_frame			= h.nx * h.ny;														// number of pixels in the image
		uint32_t n_bytes_in_binary_image 	= ceil (n_pixels_in_frame / 8.0);									// number of bytes needed to pack binary image


		// create data buffers
		uint16_t *pixvals 					= (uint16_t*)malloc (n_pixels_in_frame * sizeof(uint16_t));			// allocated max space - where all pixels are fg pixels, must calloc to init to 0
		uint8_t  *binaryImage 				= (uint8_t*) malloc (n_bytes_in_binary_image * sizeof(uint8_t));	// every bit in binaryImage has to be reset for every frame in the loop below
		uint8_t  *packedPixvals 			= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t) * 1.5);
		uint8_t  *compressedBinaryImage 	= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
		uint8_t  *compressedPackedPixvals	= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t) * 1.5);

		
		// create part file
		char* part_filename = makePartFilename(process_id, original_filename, 1);
		FILE *fp = fopen (part_filename, "wb");

	
		// serialize part header
		fwrite (&process_id, sizeof(uint8_t), 1, fp);
		
		
		uint32_t z;
		uint32_t n_bytes_in_packed_pixvals;
		uint32_t n_compressed_bytes_1, n_compressed_bytes_2;

		for (z = 0; z < h.nz; z++) {

			uint32_t frame_id = frame_start_index + z;		// absolute frame index
			
			//printf("RCT %d: 1.\n", process_id);
			//printf("RCT %d: Chunk Offset = %lu, Frame = %lu\n", process_id, frame_offset, frame_id);
			
			n_bytes_in_packed_pixvals 	= 0;
			n_compressed_bytes_1		= 0;
			n_compressed_bytes_2		= 0;
			
			reduceCompressFrame_L1 (process_id, frameBuffer, darkFrame, frame_id, z, h, input_params->rc_operation_mode, input_params->dark_threshold_epsilon,
									input_params->bit_depth, input_params->compression_scheme, input_params->compression_level, 
									pixvals, binaryImage, packedPixvals, compressedBinaryImage, compressedPackedPixvals, 
									&n_bytes_in_packed_pixvals, &n_compressed_bytes_1, &n_compressed_bytes_2, compression_time,
									0, NULL, NULL);
			//printf("RCT %d: 2.\n", process_id);
			if (input_params->rc_operation_mode 		== RC_MODE_REDUCE_ONLY) {

				fwrite (&frame_id, 					sizeof(uint32_t), 1, fp);
				fwrite (&n_bytes_in_packed_pixvals, sizeof(uint32_t), 1, fp);
				fwrite (binaryImage , 				sizeof(uint8_t),  n_bytes_in_binary_image, fp);
				fwrite (packedPixvals , 			sizeof(uint8_t),  n_bytes_in_packed_pixvals, fp);
			
			} else if (input_params->rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {

				fwrite (&frame_id, 					sizeof(uint32_t), 1, fp);
				fwrite (&n_compressed_bytes_1, 		sizeof(uint32_t), 1, fp);
				fwrite (&n_compressed_bytes_2, 		sizeof(uint32_t), 1, fp);
				fwrite (&n_bytes_in_packed_pixvals, sizeof(uint32_t), 1, fp);
				fwrite (compressedBinaryImage , 	sizeof(uint8_t),  n_compressed_bytes_1, fp);
				fwrite (compressedPackedPixvals , 	sizeof(uint8_t),  n_compressed_bytes_2, fp);
			
			}
		}


		// serialize the number of frames at the end
		fwrite (&h.nz, sizeof(uint32_t), 1, fp);

		fclose (fp);


		// clean-up
		free(pixvals);
		free(binaryImage);
		free(packedPixvals);
		free(compressedBinaryImage);
		free(compressedPackedPixvals);

		return part_filename;

}



void merge_RC1_Parts(const char* folderpath, 
					 char** 	 part_filenames, 
					 int 		 num_partfiles, 
					 RCHeader 	 *rcHeader, 
					 const char  *compressed_filename) {
	
	
	clock_t merge_start = clock();
	
	int i,j;
	
	uint32_t num_frames;
	uint32_t total_frames = 0;
	uint32_t *process_id_num_frames_map = malloc((num_partfiles)*sizeof(long));


	for (i = 0; i<num_partfiles; i++) {
		
		FILE *fp = fopen(concat(folderpath, part_filenames[i]), "rb");
		//printf("%s\n", concat(folderpath, part_filenames[i]));
		
		fseek(fp, 0, SEEK_END);
		uint32_t length = ftell(fp);
		fseek(fp, (length - 4), SEEK_SET);
		fread(&num_frames, sizeof(uint32_t), 1, fp);
		
		process_id_num_frames_map[i] = num_frames;
		total_frames += num_frames;
		
		fclose(fp);
		
		//printf("%lu\n", num_frames);
	}
	
	//printf("1\n");
	
	uint32_t *frame_process_id_map 	= malloc((total_frames)*sizeof(uint32_t));
	//uint32_t *frame_position_map 	= malloc((total_frames)*sizeof(uint32_t));
	uint32_t *frame_data_sizes_map 	= malloc((total_frames)*sizeof(uint32_t)*3);
	
	FILE *partFiles[num_partfiles];
	for (i = 0; i<num_partfiles; i++) {
		
		partFiles[i] = fopen(concat(folderpath, part_filenames[i]), "rb");
		
		//skip process id
		fseek (partFiles[i], sizeof(uint8_t), SEEK_SET);
		
		for (j = 0; j < process_id_num_frames_map[i]; j++) {
			
			uint32_t frame_id;
			fread (&frame_id, sizeof(uint32_t), 1, partFiles[i]);
			//printf("Frame ID: %d\n", frame_id);
			
			uint32_t nCompressedSize_BinaryImage;
			fread (&nCompressedSize_BinaryImage, sizeof(uint32_t), 1, partFiles[i]);
			//printf("Compressed Size 1: %d\n", nCompressedSize_BinaryImage);
			
			uint32_t nCompressedSize_Pixvals;
			fread (&nCompressedSize_Pixvals, sizeof(uint32_t), 1, partFiles[i]);
			//printf("Compressed Size 2: %d\n", nCompressedSize_Pixvals);
			
			uint32_t bytesRequiredForPacking;
			fread (&bytesRequiredForPacking, sizeof(uint32_t), 1, partFiles[i]);
			//printf("Packed Bytes: %d\n", bytesRequiredForPacking);
			
			fseek (partFiles[i], sizeof(uint8_t)*nCompressedSize_BinaryImage, SEEK_CUR);
			
			fseek (partFiles[i], sizeof(uint8_t)*nCompressedSize_Pixvals, SEEK_CUR);
			
			frame_process_id_map[frame_id] 		= i;
			//frame_position_map[frame_id] 		= j;
			frame_data_sizes_map[frame_id*3] 	= nCompressedSize_BinaryImage;
			frame_data_sizes_map[frame_id*3+1] 	= nCompressedSize_Pixvals;
			frame_data_sizes_map[frame_id*3+2] 	= bytesRequiredForPacking;
		}
		
		fclose(partFiles[i]);
	}
	
	//printf("2\n");
	
	uint32_t frame_id = 0;
	/*
	for (frame_id = 0; frame_id < total_frames; frame_id++) {
		//printf(" %" PRIu32 ", %" PRIu32 ", %" PRIu32 "\n", frame_id, frame_process_id_map[frame_id], frame_position_map[frame_id]);
		printf(" %" PRIu32 ", %" PRIu32 "\n", frame_id, frame_process_id_map[frame_id]);
	}
	*/

	compressed_filename = concat(compressed_filename, ".rc1");
	FILE *target_fp = fopen(compressed_filename, "wb");
	for (i = 0; i<num_partfiles; i++) {
		partFiles[i] = fopen(concat(folderpath, part_filenames[i]), "rb");
		// skip process_id
		fseek (partFiles[i], sizeof(uint8_t), SEEK_SET);
	}
	
	
	// write RC header to compressed_filename
	serialize_recode_header (target_fp, rcHeader);
	
	
	for (frame_id = 0; frame_id < total_frames; frame_id++) {
		// write compressed frame sizes to compressed_filename
		fwrite (&frame_data_sizes_map[frame_id*3], 	 sizeof(uint32_t), 1, target_fp);
		fwrite (&frame_data_sizes_map[frame_id*3+1], sizeof(uint32_t), 1, target_fp);
		fwrite (&frame_data_sizes_map[frame_id*3+2], sizeof(uint32_t), 1, target_fp);
	}
	
	
	// copy actual data
	uint32_t partfile_num;
	for (frame_id = 0; frame_id < total_frames; frame_id++) {
		
		partfile_num = frame_process_id_map[frame_id];
		//partfile_offset = frame_position_map[frame_id];
		
		// read a frame from partfile_num
		// skip frame_id, nCompressedSize_BinaryImage, nCompressedSize_Pixvals, bytesRequiredForPacking
		fseek (partFiles[partfile_num], sizeof(uint32_t)*4, SEEK_CUR);
			
		uint8_t *compressedBinaryImage = (uint8_t*)calloc(frame_data_sizes_map[frame_id*3], sizeof(uint8_t));
		fread (compressedBinaryImage, sizeof(uint8_t), frame_data_sizes_map[frame_id*3], partFiles[partfile_num]);
			
		uint8_t *compressedPixvals = (uint8_t*)calloc(frame_data_sizes_map[frame_id*3+1], sizeof(uint8_t));
		fread (compressedPixvals, sizeof(uint8_t), frame_data_sizes_map[frame_id*3+1], partFiles[partfile_num]);
		
		// write frame to compressed_filename
		fwrite (compressedBinaryImage , sizeof(uint8_t), frame_data_sizes_map[frame_id*3], target_fp);
		fwrite (compressedPixvals , sizeof(uint8_t), frame_data_sizes_map[frame_id*3+1], target_fp);
	}
	
	for (i = 0; i<num_partfiles; i++) {
		fclose(partFiles[i]);
	}
	fclose(target_fp);
	
	free(process_id_num_frames_map);
	free(frame_process_id_map);
	//free(frame_position_map);
	free(frame_data_sizes_map);
	
	printf("111\n");
	
	clock_t merge_end = clock();
	float merge_time = (merge_end - merge_start) * 1000.0 / CLOCKS_PER_SEC;
	printf("Merge Time: %f\n", merge_time);
}


void decompressExpand_L1 (	const char* compressed_filename, 
							uint16_t	**frameBuffer, 
							RCHeader 	**header) {
	
	
	FILE* fp = fopen(compressed_filename,"rb");
	
	// read RC header from compressed file
	parse_recode_header (fp, header);
	
	// read compressed data lengths
	uint32_t frame_id = 0;
	uint32_t *frame_data_sizes_map 	= malloc(((*header)->nz)*sizeof(uint32_t)*3);
	for (frame_id = 0; frame_id < (*header)->nz; frame_id++) {
		fread (&frame_data_sizes_map[frame_id*3], 	sizeof(uint32_t), 1, fp);
		fread (&frame_data_sizes_map[frame_id*3+1], sizeof(uint32_t), 1, fp);
		fread (&frame_data_sizes_map[frame_id*3+2], sizeof(uint32_t), 1, fp);
	}
	
	/*
	for (frame_id = 0; frame_id < (*header)->nz; frame_id++) {
		printf( " %" PRIu32 ": %" PRIu32 ", %" PRIu32 ", %" PRIu32 "\n", 
				frame_id, frame_data_sizes_map[frame_id*3], 
				frame_data_sizes_map[frame_id*3+1], 
				frame_data_sizes_map[frame_id*3+2]);
	}
	*/
	
	uint32_t n_pixels_in_frame			= (*header)->nx * (*header)->ny;
	uint32_t n_bytes_in_binary_image 	= ceil((n_pixels_in_frame*1.0) / 8.0);
	
	
	printf("Allocating memory.\n");
	uint8_t bytes_per_pixel = ceil( ((*header)->bit_depth*1.0) / 8.0 );
	uint8_t bits_per_pixel 	= bytes_per_pixel * 8;
	*frameBuffer = (uint16_t*)  malloc ((*header)->nx * (*header)->ny * (*header)->nz * bytes_per_pixel);
	
	uint16_t extracted_pixval = 0;
	
	uint64_t pow2_lookup[64];
	uint8_t	 t;
	for (t = 0; t < 64; t++) {
		pow2_lookup[t] = pow(2,t);
		//printf ("Power 2 of %d = %"PRIu64"\n", t, pow2_lookup[t]);
	}
	
	// read frames
	for (frame_id = 0; frame_id < (*header)->nz; frame_id++) {
		
		//printf("Frame ID: %d\n", frame_id);
		
		uint32_t n_compressed_bytes_binary_img 	= frame_data_sizes_map[frame_id*3];
		uint32_t n_compressed_bytes_pixvals 	= frame_data_sizes_map[frame_id*3+1];
		uint32_t n_bytes_in_packed_pixvals 		= frame_data_sizes_map[frame_id*3+2];
		
		uint8_t *compressedBinaryImage 	= (uint8_t*) calloc (frame_data_sizes_map[frame_id*3], sizeof(uint8_t));
		uint8_t *compressedPixvals 		= (uint8_t*) calloc (n_compressed_bytes_pixvals, 		   sizeof(uint8_t));
		
		fread (compressedBinaryImage, sizeof(uint8_t), n_compressed_bytes_binary_img, fp);
		fread (compressedPixvals, 	  sizeof(uint8_t), n_compressed_bytes_pixvals, 	fp);
		
		uint8_t* deCompressedBinaryImage = (uint8_t*) calloc (n_bytes_in_binary_image, 	 sizeof(uint8_t));
		uint8_t* deCompressedPixvals 	 = (uint8_t*) calloc (n_bytes_in_packed_pixvals, sizeof(uint8_t));
		
		decompress_stream (-1, -1, compressedBinaryImage, deCompressedBinaryImage, n_compressed_bytes_binary_img, n_bytes_in_binary_image);
		decompress_stream (-1, -1, compressedPixvals, 	  deCompressedPixvals, 	   n_compressed_bytes_pixvals, 	  n_bytes_in_packed_pixvals);
		
		//printf("Ready... Set... Go.\n");
		uint32_t frame_start_linear_index = n_pixels_in_frame * frame_id;
		//printf ("frame_start_linear_index = %lu\n", frame_start_linear_index);
		
		uint8_t n, t, nth_bit_of_pth_pixval;
		uint32_t row, col, linear_pixel_index, pixel_bit_index_frame, pixel_bit_index_dataset;
		uint32_t fg_pixel_count = 0;
		
		uint32_t nx 		= (*header)->nx;
		uint32_t ny 		= (*header)->ny;
		uint32_t bit_depth 	= (*header)->bit_depth;
		
		//printf("Ready... Set... Go.");
		uint32_t n_fg_pixels = 0;
		for (row = 0; row < ny; row++) {
			for (col = 0; col < nx; col++) {
				linear_pixel_index = row * nx + col;
				
				/*
				if ( CheckBit(deCompressedBinaryImage, linear_pixel_index) > 0 ) {
					
					(*frameBuffer)[frame_start_linear_index + linear_pixel_index] = 1;
					//printf("row: %hu, col: %hu\n", row, col);
					n_fg_pixels++;
					
				} else {
					
					(*frameBuffer)[frame_start_linear_index + linear_pixel_index] = 0;
					
				}
				*/
				
				//printf ("linear_pixel_index = %lu\n", linear_pixel_index);
				if ( CheckBit(deCompressedBinaryImage, linear_pixel_index) > 0 ) {
					
					extracted_pixval = 0;
					
					// Assumes LITTLE-ENDIAN Byte Order of pixval
					for (n=0; n<bit_depth; n++) {
						
						pixel_bit_index_frame 	= fg_pixel_count*bit_depth + n;
						pixel_bit_index_dataset = (frame_start_linear_index + linear_pixel_index)*bits_per_pixel + n;
						
						//printf ("pixel_bit_index_frame = %lu\n", pixel_bit_index_frame);
						//printf ("pixel_bit_index_dataset = %lu\n", pixel_bit_index_dataset);
						
						if (CheckBit(deCompressedPixvals, pixel_bit_index_frame) > 0) {
							extracted_pixval += pow2_lookup[n];
						}
					}
					
					(*frameBuffer)[frame_start_linear_index + linear_pixel_index] = extracted_pixval;
					
					//printf(" ");
					//printf("count: %hu, row: %hu, col: %hu, index = %"PRIu32" val = %hu\n", 
					//		fg_pixel_count, row, col, linear_pixel_index, extracted_pixval);

					fg_pixel_count++;
					
				} else {
					
					(*frameBuffer)[frame_start_linear_index + linear_pixel_index] = 0;
					
				}
				
			}
		}
		
		/*
		uint32_t cnt = 0;
		//uint32_t row, col, linear_pixel_index;
		for (row = 0; row < ny; row++) {
			for (col = 0; col < nx; col++) {
				linear_pixel_index = row * nx + col;
				if ( CheckBit(deCompressedBinaryImage, linear_pixel_index) > 0 ) {
					printf("count: %hu, row: %hu, col: %hu, index: %lu, val = %hu\n", 
							cnt, row, col, linear_pixel_index, (*frameBuffer)[linear_pixel_index]);
					cnt++;
				}
			}
		}
		*/
	
		//printf("bytesRequiredForPacking: %d\n", bytesRequiredForPacking);
		printf ("n_fg_pixels = %lu\n", fg_pixel_count);
		
	}
	printf(" Done.\n");
}

