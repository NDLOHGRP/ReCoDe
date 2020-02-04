


// Used in L1
float binarize_and_get_pixvals (uint8_t	 process_id,
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
	
	//recode_print("RCT %d: No. of foreground pixels = %lu\n", process_id, *n_fg_pixels);

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}

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

	uint64_t sz_packedPixval = (int)ceil(n_fg_pixels*(bit_depth / 8.0));
	for (p = 0; p < sz_packedPixval; p++) {
		packedPixvals[p] = 0;
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
			linear_index = p*bit_depth + n;
			if (nth_bit_of_pth_pixval == 1) {
				SetBit(packedPixvals, linear_index);
			}
			// setting packedPixvals to 0 in the for loop (above) is faster than using ClearBit
			/*
			else {
				ClearBit(packedPixvals, linear_index);
			}
			*/
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
		
	float thresh_time 			= binarize_and_get_pixvals (process_id, frameBuffer + frame_start_index, darkFrame, epsilon_s, h, 
															binaryImage, pixvals, &n_fg_pixels, &data_min, &data_max);
	*n_bytes_in_packed_pixvals 	= ceil((n_fg_pixels * (bit_depth))/8);
	float packing_time 			= scale_and_pack_pixvals (process_id, pixvals, n_fg_pixels, data_min, data_max, bit_depth, packedPixvals);
	//printf("packing_time=%f\n", packing_time);
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
	//printf("RCT %d: 1.2.1\n", process_id);
	float 	 compression_time_1 = compress_stream (compression_scheme, compression_level, binaryImage,   n_bytes_in_binary_image, n_compressed_bytes_1, compressedBinaryImage);
	//printf("RCT %d: 1.2.2\n", process_id);
	float 	 compression_time_2 = compress_stream (compression_scheme, compression_level, packedPixvals, *n_bytes_in_packed_pixvals, n_compressed_bytes_2, compressedPackedPixvals);
	
	float 	 compression_time 	= compression_time_1 + compression_time_2;
	
	return compression_time;
}

float decompressFrame_L1 (	uint8_t	 process_id,
							DataSize h, 
							uint8_t  compression_scheme,
							uint8_t  compression_level,
							uint8_t	 *binaryImage,
							uint8_t  *packedPixvals,
							uint32_t *n_bytes_in_binary_image,
							uint32_t *n_bytes_in_packed_pixvals,
							uint8_t  *compressedBinaryImage,
							uint8_t  *compressedPackedPixvals,
							uint32_t *n_compressed_bytes_binary_img,
							uint32_t *n_compressed_bytes_packed_pixvals
							) {
	
	/*
	========================================================================================== 
	Decompress
	========================================================================================== 
	*/
	uint8_t  *temp_1			 		= (uint8_t*)malloc (*n_bytes_in_binary_image);
	uint8_t  *temp_2			 		= (uint8_t*)malloc (*n_bytes_in_packed_pixvals);

	float 	 decompression_time_1 		= decompress_stream (compression_scheme, compression_level, 
															 compressedBinaryImage, temp_1, 
															 *n_compressed_bytes_binary_img, *n_bytes_in_binary_image);
	printf ("decompressed binary image.\n");
	
	float 	 decompression_time_2 		= decompress_stream (compression_scheme, compression_level, 
															 compressedPackedPixvals, temp_2, 
															 *n_compressed_bytes_packed_pixvals, *n_bytes_in_packed_pixvals);
	printf ("decompressed pixvals.\n");
	
	float 	 decompression_time 		= decompression_time_1 + decompression_time_2;
	
	
	/*
	========================================================================================== 
	Validate Decompressed Data
	========================================================================================== 
	*/
	uint32_t i;
	float percent_mismatch = 0.0;
	for (i = 0; i < *n_bytes_in_binary_image; i++) {
		if (temp_1[i] != binaryImage[i]) {
			percent_mismatch++;
		}
	}
	if (percent_mismatch > -1) {
		percent_mismatch /= *n_bytes_in_binary_image;
		printf("Compression Scheme %d: Mismatch % (binary image) = %.2f\n", compression_scheme, percent_mismatch);
	}

	percent_mismatch = 0.0;
	for (i = 0; i < *n_bytes_in_packed_pixvals; i++) {
		if (temp_2[i] != packedPixvals[i]) {
			percent_mismatch++;
		}
	}
	if (percent_mismatch > -1) {
		percent_mismatch /= *n_bytes_in_packed_pixvals;
		printf("Compression Scheme %d: Mismatch % (packed pixvals) = %.2f\n", compression_scheme, percent_mismatch);
	}
	
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

		uint32_t n_bytes_in_binary_image = ceil((h.nx * h.ny) / 8.0);
								
		//printf("RCT %d: compressed_frame[0] = %d\n", process_id, compressed_frame[0]);
								
		//printf("RCT %d: 1.1.\n", process_id);
		float reduction_time  	= reduceFrame_L1 (process_id, frameBuffer, darkFrame, z, h, 
												  epsilon_s, bit_depth, 
												  pixvals, binaryImage, packedPixvals, n_bytes_in_packed_pixvals);

		float compression_time  	= 0.0;
		float decompression_time  	= 0.0;

		if (rc_operation_mode == RC_MODE_REDUCE_COMPRESS) {
			//printf("RCT %d: 1.2.\n", process_id);
			compression_time  	= compressFrame_L1 (process_id, binaryImage, packedPixvals, n_bytes_in_packed_pixvals, h,
													compression_scheme, compression_level,
						  							compressedBinaryImage, compressedPackedPixvals, 
						  							n_compressed_bytes_1, n_compressed_bytes_2);
			//printf("RCT %d: 1.3.\n", process_id);

			/*
			decompression_time	= decompressFrame_L1 (	process_id, h,
														compression_scheme, compression_level,
														binaryImage, packedPixvals, 
														&n_bytes_in_binary_image, n_bytes_in_packed_pixvals,
						  								compressedBinaryImage, compressedPackedPixvals, 
						  								n_compressed_bytes_1, n_compressed_bytes_2);
			*/
		}
		

		//recode_print ("RCT %d: Frame ID = %lu, n_bytes_in_packed_pixvals = %lu\n", process_id, frame_id, *n_bytes_in_packed_pixvals);

		run_metrics[0] += reduction_time;
		run_metrics[1] += compression_time;
		//run_metrics[2] += (n_bytes_in_binary_image + *n_bytes_in_packed_pixvals);
        run_metrics[2] += *n_bytes_in_packed_pixvals;                                   // for unreduced compression test
		//run_metrics[3] += (*n_compressed_bytes_1 + *n_compressed_bytes_2);
		run_metrics[3] += (*n_compressed_bytes_2);                                      // for unreduced compression test
		run_metrics[4] += decompression_time;
		
		
		/*
		========================================================================================== 
		Copy compressed data into return buffer if requested. Used in on-the-fly compression.
		========================================================================================== 
		*/
		if (copy_compressed_frame_to_return_buffer == 1) {
			
			if (rc_operation_mode 			== RC_MODE_REDUCE_ONLY) {
			
				*compressed_frame_length 	= sizeof(uint8_t) * (n_bytes_in_binary_image) + sizeof(uint8_t) * (*n_bytes_in_packed_pixvals) + 2 * sizeof(uint32_t);

				memcpy(compressed_frame,	&frame_id, 					sizeof(uint32_t));
				memcpy(compressed_frame+4,	n_bytes_in_packed_pixvals, 	sizeof(uint32_t));
				memcpy(compressed_frame+8,	binaryImage, 				sizeof(uint8_t) * n_bytes_in_binary_image);
				memcpy(compressed_frame+8+n_bytes_in_binary_image, packedPixvals, sizeof(uint8_t) * (*n_bytes_in_packed_pixvals));
				
				///printf("Copying to return buffer in Reduce Only Mode.\n");
			
			} else if (rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {
				
				*compressed_frame_length 	= sizeof(uint8_t) * (*n_compressed_bytes_1) + sizeof(uint8_t) * (*n_compressed_bytes_2) + 4 * sizeof(uint32_t);
				
				memcpy(compressed_frame,		&frame_id, 					sizeof(uint32_t));
				memcpy(compressed_frame+4,		n_compressed_bytes_1, 		sizeof(uint32_t));
				memcpy(compressed_frame+8,		n_compressed_bytes_2, 		sizeof(uint32_t));
				memcpy(compressed_frame+12,		n_bytes_in_packed_pixvals, 	sizeof(uint32_t));
				memcpy(compressed_frame+16,		compressedBinaryImage, 		sizeof(uint8_t)*(*n_compressed_bytes_1));
				memcpy(compressed_frame+16+(*n_compressed_bytes_1), compressedPackedPixvals, 	sizeof(uint8_t)*(*n_compressed_bytes_2));
				
				///printf("Copying to return buffer in Reduce Compress Mode.\n");
			}
		}
		
		return;
	
}



char* reduceCompress_L1 (uint8_t	 process_id, 
						 const char  *original_filename, 
						 const char  *out_foldername, 
						 uint16_t    *frameBuffer, 				// points to the start of this thread's chunk
						 uint16_t    *darkFrame, 
						 uint32_t    frame_start_index, 
						 DataSize    h, 
						 InputParams *input_params,
						 RCHeader 	 *rcHeader,
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
		char *part_filename = makePartFilename(process_id, original_filename, 1);
		FILE *fp = fopen (concat(out_foldername, part_filename), "wb");

		// serialize header
		serialize_recode_header(fp, rcHeader);
		
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
		// fwrite (&h.nz, sizeof(uint32_t), 1, fp);

		// serialize the number of frames in the header
		fseek(fp, 17, SEEK_SET);
		fwrite(&h.nz, sizeof(uint32_t), 1, fp);
		fseek(fp, 277, SEEK_SET);
		fwrite(&process_id, sizeof(uint8_t), 1, fp);
		fclose(fp);

		// clean-up
		free(pixvals);
		free(binaryImage);
		free(packedPixvals);
		free(compressedBinaryImage);
		free(compressedPackedPixvals);

		return part_filename;

}



char* merge_RC1_Parts (	const char   *folderpath, 
						char         **part_filenames, 
						InputParams  *input_params
) {
	
	const int num_partfiles = input_params->num_threads;
	clock_t merge_start = clock();
	
	int i,j;
	
	uint32_t total_frames = 0;
	uint32_t *process_id_num_frames_map = (uint32_t*)malloc((num_partfiles)*sizeof(long));

	RCHeader *header = (RCHeader *)malloc(sizeof(RCHeader));
	for (i = 0; i<num_partfiles; i++) {
		FILE *fp = fopen(concat(folderpath, part_filenames[i]), "rb");
		recode_print("%s\n", concat(folderpath, part_filenames[i]));
		parse_recode_header(fp, &header);
		process_id_num_frames_map[i] = header->nz;
		total_frames += header->nz;
		fclose(fp);
		print_recode_header(header);
		recode_print("%lu\n", total_frames);
	}
	
	//printf("1\n");
	
	int32_t *frame_process_id_map 	= (int32_t*)malloc((total_frames)*sizeof(int32_t));
	uint32_t frame_id = 0;
	for (frame_id = 0; frame_id < total_frames; frame_id++) {
		frame_process_id_map[frame_id] = -1;
	}
	//uint32_t *frame_data_sizes_map 	= (uint32_t*)malloc((total_frames)*sizeof(uint32_t)*3);
	uint32_t *frame_data_sizes_map = (uint32_t*)calloc((total_frames) * 3, sizeof(uint32_t));

	FILE **partFiles = (FILE**)malloc((num_partfiles)*sizeof(FILE*));
	for (i = 0; i<num_partfiles; i++) {
		
		partFiles[i] = fopen(concat(folderpath, part_filenames[i]), "rb");
		
		//skip header
		fseek(partFiles[i], RC_HEADER_LENGTH, SEEK_SET);
		
		for (j = 0; j < process_id_num_frames_map[i]; j++) {
			
			uint32_t frame_id;
			fread (&frame_id, sizeof(uint32_t), 1, partFiles[i]);
			
			uint32_t nCompressedSize_BinaryImage;
			fread (&nCompressedSize_BinaryImage, sizeof(uint32_t), 1, partFiles[i]);
			
			uint32_t nCompressedSize_Pixvals;
			fread (&nCompressedSize_Pixvals, sizeof(uint32_t), 1, partFiles[i]);
			
			uint32_t bytesRequiredForPacking;
			fread (&bytesRequiredForPacking, sizeof(uint32_t), 1, partFiles[i]);
			
			fseek (partFiles[i], sizeof(uint8_t)*nCompressedSize_BinaryImage, SEEK_CUR);
			fseek (partFiles[i], sizeof(uint8_t)*nCompressedSize_Pixvals, SEEK_CUR);
			
			frame_process_id_map[frame_id] 		= i;
			frame_data_sizes_map[frame_id*3] 	= nCompressedSize_BinaryImage;
			frame_data_sizes_map[frame_id*3+1] 	= nCompressedSize_Pixvals;
			frame_data_sizes_map[frame_id*3+2] 	= bytesRequiredForPacking;

			/*
			printf("Frame ID: %d\n", frame_id);
			printf("Compressed Size 1: %d\n", nCompressedSize_BinaryImage);
			printf("Compressed Size 2: %d\n", nCompressedSize_Pixvals);
			printf("Packed Bytes: %d\n", bytesRequiredForPacking);
			*/
		}
		
		fclose(partFiles[i]);
	}
	
	const char *compressed_filename = concat(folderpath, concat((const char*)&(header->source_file_name), ".rc1"));
	FILE *target_fp = fopen(compressed_filename, "wb");
	for (i = 0; i<num_partfiles; i++) {
		partFiles[i] = fopen(concat(folderpath, part_filenames[i]), "rb");
		//skip header
		fseek(partFiles[i], RC_HEADER_LENGTH, SEEK_SET);
	}
	
	// write header to compressed_filename
	// serialize the number of frames in the header
	header->nz = total_frames;
	serialize_recode_header(target_fp, header);
	
	for (frame_id = 0; frame_id < total_frames; frame_id++) {
		// write compressed frame sizes to compressed_filename
		fwrite (&frame_data_sizes_map[frame_id*3], 	 sizeof(uint32_t), 1, target_fp);
		fwrite (&frame_data_sizes_map[frame_id*3+1], sizeof(uint32_t), 1, target_fp);
		fwrite (&frame_data_sizes_map[frame_id*3+2], sizeof(uint32_t), 1, target_fp);
	}

	// copy actual data
	uint32_t partfile_num;
	uint32_t n_bytes_in_image = ceil(input_params->num_rows * input_params->num_cols * input_params->bit_depth / 8.0);
	uint32_t n_bytes_in_binary_image = ceil(input_params->num_rows * input_params->num_cols / 8.0);
	uint8_t *compressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	for (frame_id = 0; frame_id < total_frames; frame_id++) {
		
		// get the part file that contains the frame
		partfile_num = frame_process_id_map[frame_id];
		if (partfile_num == -1) {
			printf("Missing Frame: %d\n", frame_id);
			continue;
		}
		
		// skip frame_id, nCompressedSize_BinaryImage, nCompressedSize_Pixvals, bytesRequiredForPacking
		fseek (partFiles[partfile_num], sizeof(uint32_t) * 4, SEEK_CUR);
			
		// read frame data from part file
		fread (compressedBinaryImage, sizeof(uint8_t), frame_data_sizes_map[frame_id*3], partFiles[partfile_num]);
		fread (compressedPixvals, sizeof(uint8_t), frame_data_sizes_map[frame_id*3+1], partFiles[partfile_num]);
		
		// write frame data to compressed_filename
		fwrite (compressedBinaryImage , sizeof(uint8_t), frame_data_sizes_map[frame_id*3], target_fp);
		fwrite (compressedPixvals , sizeof(uint8_t), frame_data_sizes_map[frame_id*3+1], target_fp);

	}
	
	for (i = 0; i<num_partfiles; i++) {
		fclose(partFiles[i]);
	}
	fclose(target_fp);

	free(header);
	free(compressedBinaryImage);
	free(compressedPixvals);
	free(process_id_num_frames_map);
	free(frame_process_id_map);
	free(frame_data_sizes_map);
	free(partFiles);

	clock_t merge_end = clock();
	float merge_time = (merge_end - merge_start) * 1000.0 / CLOCKS_PER_SEC;
	printf("Merge Time: %f\n", merge_time);
	
	return (char*)compressed_filename;
}

/*
This function is provided for testing purposes only. No external pythonic calls to this function are available. Decompression to only sparse format is supported for external calls.
*/
void decompressExpand_L1_Reduced_Compressed ( FILE* rc_fp, const char *out_fname, RCHeader *header) {

	// fp is assumed to point to end of header, as header was read using the same fp in recode.cpp
	
	uint32_t nx = header->nx;
	uint32_t ny = header->ny;
	uint32_t nz = header->nz;
	uint8_t byte_depth = ceil((header->bit_depth*1.0) / 8.0);

	uint32_t n_pixels_in_frame = nx * ny;														// number of pixels in the image
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);							// number of bytes needed to pack binary image
	uint32_t n_bytes_in_image = ceil(n_pixels_in_frame * byte_depth);							// number of bytes needed to hold pixvals

	uint64_t sz_frameBuffer = n_pixels_in_frame * byte_depth;
	uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);

	uint32_t frame_id = 0;
	uint32_t *n_compressed_bytes_in_binary_image = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_compressed_bytes_in_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_bytes_in_packed_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	for (frame_id = 0; frame_id < header->nz; frame_id++) {
		fread( &n_compressed_bytes_in_binary_image [frame_id], sizeof(uint32_t), 1, rc_fp);
		fread( &n_compressed_bytes_in_pixvals [frame_id], sizeof(uint32_t), 1, rc_fp);
		fread( &n_bytes_in_packed_pixvals [frame_id], sizeof(uint32_t), 1, rc_fp);
	}

	uint8_t *compressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *deCompressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));

	uint8_t *compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	uint8_t *deCompressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));

	uint64_t pow2_lookup_table[64];
	uint8_t	 t;
	for (t = 0; t < 64; t++) {
		pow2_lookup_table[t] = pow(2, t);
	}

	// read frames
	uint8_t n;
	uint16_t extracted_pixval;

	FILE *out_fp = fopen(out_fname, "wb+");
	serialize_recode_header(out_fp, header);

	for (frame_id = 0; frame_id < nz; frame_id++) {
		fread(compressedBinaryImage, sizeof(uint8_t), n_compressed_bytes_in_binary_image[frame_id], rc_fp);
		fread(compressedPixvals, sizeof(uint8_t), n_compressed_bytes_in_pixvals[frame_id], rc_fp);

		decompress_stream(-1, -1, compressedBinaryImage, deCompressedBinaryImage, n_compressed_bytes_in_binary_image[frame_id], n_bytes_in_binary_image);
		decompress_stream(-1, -1, compressedPixvals, deCompressedPixvals, n_compressed_bytes_in_pixvals[frame_id], n_bytes_in_packed_pixvals[frame_id]);

		uint32_t frame_start_linear_index = n_pixels_in_frame * frame_id;
		uint32_t row, col, linear_pixel_index, pixel_bit_index_frame;
		//uint64_t pixel_bit_index_dataset;
		uint64_t n_fg_pixels = 0;
		for (row = 0; row < ny; row++) {
			for (col = 0; col < nx; col++) {
				linear_pixel_index = row * nx + col;
				if (CheckBit(deCompressedBinaryImage, linear_pixel_index) > 0) {
					// unpack pixval
					extracted_pixval = 0;
					// Assumes LITTLE-ENDIAN Byte Order of pixval
					for (n = 0; n<header->bit_depth; n++) {
						pixel_bit_index_frame = n_fg_pixels*header->bit_depth + n;
						//pixel_bit_index_dataset = (frame_start_linear_index + linear_pixel_index)*(byte_depth * 8) + n;
						if (CheckBit(deCompressedPixvals, pixel_bit_index_frame) > 0) {
							extracted_pixval += pow2_lookup_table[n];
						}
					}
					frameBuffer[linear_pixel_index] = extracted_pixval;
					n_fg_pixels++;
				}
				else {
					frameBuffer[linear_pixel_index] = 0;
				}
			}
		}
		fwrite(frameBuffer, sizeof(uint16_t), n_pixels_in_frame, out_fp);
		printf("Decoded Frame %" PRIu32 " with %" PRIu64 " foreground pixels\n", frame_id, n_fg_pixels);
	}
	
	fclose(out_fp);
	recode_print("Done.\n");
}


void decompressExpand_L1_Reduced_Compressed_Sparse (FILE* rc_fp, const char *out_fname, RCHeader *header) {

	// fp is assumed to point to end of header, as header was read using the same fp in recode.cpp

	uint32_t nx = header->nx;
	uint32_t ny = header->ny;
	uint32_t nz = header->nz;
	uint8_t byte_depth = ceil((header->bit_depth*1.0) / 8.0);

	uint32_t n_pixels_in_frame = nx * ny;														// number of pixels in the image
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);							// number of bytes needed to pack binary image
	uint32_t n_bytes_in_image = ceil(n_pixels_in_frame * byte_depth);							// number of bytes needed to hold pixvals

	uint64_t sz_frameBuffer = n_pixels_in_frame * 3;
	uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);

	uint32_t frame_id = 0;
	uint32_t *n_compressed_bytes_in_binary_image = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_compressed_bytes_in_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_bytes_in_packed_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	for (frame_id = 0; frame_id < header->nz; frame_id++) {
		fread(&n_compressed_bytes_in_binary_image[frame_id], sizeof(uint32_t), 1, rc_fp);
		fread(&n_compressed_bytes_in_pixvals[frame_id], sizeof(uint32_t), 1, rc_fp);
		fread(&n_bytes_in_packed_pixvals[frame_id], sizeof(uint32_t), 1, rc_fp);
	}

	uint8_t *compressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *deCompressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));

	uint8_t *compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	uint8_t *deCompressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));

	uint64_t pow2_lookup_table[64];
	uint8_t	 t;
	for (t = 0; t < 64; t++) {
		pow2_lookup_table[t] = pow(2, t);
	}

	// read frames
	uint8_t n;
	uint16_t extracted_pixval;

	FILE *out_fp = fopen(out_fname, "wb+");
	serialize_recode_header(out_fp, header);

	for (frame_id = 0; frame_id < nz; frame_id++) {
		fread(compressedBinaryImage, sizeof(uint8_t), n_compressed_bytes_in_binary_image[frame_id], rc_fp);
		fread(compressedPixvals, sizeof(uint8_t), n_compressed_bytes_in_pixvals[frame_id], rc_fp);

		decompress_stream(-1, -1, compressedBinaryImage, deCompressedBinaryImage, n_compressed_bytes_in_binary_image[frame_id], n_bytes_in_binary_image);
		decompress_stream(-1, -1, compressedPixvals, deCompressedPixvals, n_compressed_bytes_in_pixvals[frame_id], n_bytes_in_packed_pixvals[frame_id]);

		//uint32_t frame_start_linear_index = n_pixels_in_frame * frame_id;
		uint16_t row, col;
		uint32_t linear_pixel_index, pixel_bit_index_frame;
		//uint64_t pixel_bit_index_dataset;
		uint64_t n_fg_pixels = 0;
		for (row = 0; row < ny; row++) {
			for (col = 0; col < nx; col++) {
				linear_pixel_index = row * nx + col;
				if (CheckBit(deCompressedBinaryImage, linear_pixel_index) > 0) {
					// unpack pixval
					extracted_pixval = 0;
					// Assumes LITTLE-ENDIAN Byte Order of pixval
					for (n = 0; n<header->bit_depth; n++) {
						pixel_bit_index_frame = n_fg_pixels*header->bit_depth + n;
						//pixel_bit_index_dataset = (frame_start_linear_index + linear_pixel_index)*(byte_depth * 8) + n;
						if (CheckBit(deCompressedPixvals, pixel_bit_index_frame) > 0) {
							extracted_pixval += pow2_lookup_table[n];
						}
					}
					frameBuffer[n_fg_pixels * 3] = row;
					frameBuffer[n_fg_pixels * 3 + 1] = col;
					frameBuffer[n_fg_pixels * 3 + 2] = extracted_pixval;
					n_fg_pixels++;
				}
			}
		}
		fwrite(&n_fg_pixels, sizeof(uint32_t), 1, out_fp);
		fwrite(frameBuffer, sizeof(uint16_t), n_fg_pixels * 3, out_fp);
		printf("Decoded Frame %" PRIu32 " with %" PRIu64 " foreground pixels\n", frame_id, n_fg_pixels);
	}

	fclose(out_fp);
	recode_print("Done.\n");
}



/*
this function is called by pyrecode_c.py
*/
int64_t decompressExpand_L1_Reduced_Compressed_Frame_Sparse(
	FILE* rc_fp, uint16_t nx, uint16_t ny, uint8_t bit_depth, 
	uint32_t n_compressed_bytes_in_binary_image, uint32_t n_compressed_bytes_in_pixvals, uint32_t n_bytes_in_packed_pixvals, uint32_t n_bytes_in_binary_image,
	uint8_t *compressedBinaryImage, uint8_t *deCompressedBinaryImage, uint8_t *compressedPixvals, uint8_t *deCompressedPixvals,
	uint64_t *pow2_lookup_table,
	uint16_t *frameBuffer,
	uint8_t is_intermediate_file
) {

	/*
	rc_fp is assumed to point to beginning of a frame. Appropriate seek is handled by callers: pyrecode_c.py
	*/
	size_t t;
	uint32_t frame_id;
	if (is_intermediate_file == 1) {
		// skip frame_id, nCompressedSize_BinaryImage, nCompressedSize_Pixvals, bytesRequiredForPacking
		t = fread(&frame_id, sizeof(uint32_t), 1, rc_fp);
		if (t != 1) {
			return -1;
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
	printf("L1.h: Processing Frame % PRIu32 \n", frame_id);
	printf("nx = %" PRIu32 "\n", nx);
	printf("ny = %" PRIu32 "\n", ny);
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
				//printf("Row = %" PRIu16 ", Col = %" PRIu16 ", Value = %" PRIu16 ", foreground pixel = %" PRIu64 "\n", row, col, extracted_pixval, n_fg_pixels);
				n_fg_pixels++;
			}
		}
	}
	
	printf("Decoded Frame with %" PRIu64 " foreground pixels\n", n_fg_pixels);
	return (int64_t)n_fg_pixels;
}


/*
handle for testing decompressExpand_L1_Reduced_Compressed_Frame_Sparse
*/
void _h_decompressExpand_L1_Reduced_Compressed_Sparse(FILE* rc_fp, const char *out_fname, RCHeader *header) {

	uint16_t nx = header->nx;
	uint16_t ny = header->ny;
	uint32_t nz = header->nz;
	uint8_t byte_depth = ceil((header->bit_depth*1.0) / 8.0);

	uint32_t n_pixels_in_frame = (uint32_t)nx * (uint32_t)ny;					// number of pixels in the image
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);			// number of bytes needed to pack binary image
	uint32_t n_bytes_in_image = ceil(n_pixels_in_frame * byte_depth);			// number of bytes needed to hold pixvals

	uint64_t sz_frameBuffer = n_pixels_in_frame * 3;
	uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);

	uint32_t frame_id = 0;
	uint32_t *n_compressed_bytes_in_binary_image = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_compressed_bytes_in_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_bytes_in_packed_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	for (frame_id = 0; frame_id < header->nz; frame_id++) {
		fread(&n_compressed_bytes_in_binary_image[frame_id], sizeof(uint32_t), 1, rc_fp);
		fread(&n_compressed_bytes_in_pixvals[frame_id], sizeof(uint32_t), 1, rc_fp);
		fread(&n_bytes_in_packed_pixvals[frame_id], sizeof(uint32_t), 1, rc_fp);
	}

	uint8_t *compressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *deCompressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));

	uint8_t *compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	uint8_t *deCompressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));

	uint64_t pow2_lookup_table[64];
	uint8_t	 t;
	for (t = 0; t < 64; t++) {
		pow2_lookup_table[t] = pow(2, t);
	}

	for (int i = 0; i < 3; i++) {
		decompressExpand_L1_Reduced_Compressed_Frame_Sparse(
			rc_fp, nx, ny, header->bit_depth,
			n_compressed_bytes_in_binary_image[i], n_compressed_bytes_in_pixvals[i], n_bytes_in_packed_pixvals[i], n_bytes_in_binary_image,
			compressedBinaryImage, deCompressedBinaryImage, compressedPixvals, deCompressedPixvals,
			pow2_lookup_table,
			frameBuffer,
			0
		);
	}
}

/*
handle for testing decompressExpand_L1_Reduced_Compressed_Frame_Sparse on intermediate files
*/
void _h_decompressExpand_L1_Reduced_Compressed_Sparse_Intermediate (FILE* rc_fp, const char *out_fname, RCHeader *header) {

	uint16_t nx = header->nx;
	uint16_t ny = header->ny;
	uint32_t nz = header->nz;
	uint8_t byte_depth = ceil((header->bit_depth*1.0) / 8.0);

	uint32_t n_pixels_in_frame = (uint32_t)nx * (uint32_t)ny;					// number of pixels in the image
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);			// number of bytes needed to pack binary image
	uint32_t n_bytes_in_image = ceil(n_pixels_in_frame * byte_depth);			// number of bytes needed to hold pixvals

	uint64_t sz_frameBuffer = n_pixels_in_frame * 3;
	uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);

	uint32_t frame_id = 0;
	uint32_t *n_compressed_bytes_in_binary_image = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_compressed_bytes_in_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));
	uint32_t *n_bytes_in_packed_pixvals = (uint32_t *)malloc((header->nz) * sizeof(uint32_t));

	uint8_t *compressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));
	uint8_t *deCompressedBinaryImage = (uint8_t*)calloc(n_bytes_in_binary_image, sizeof(uint8_t));

	uint8_t *compressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));
	uint8_t *deCompressedPixvals = (uint8_t*)calloc(n_bytes_in_image, sizeof(uint8_t));

	uint64_t pow2_lookup_table[64];
	uint8_t	 t;
	for (t = 0; t < 64; t++) {
		pow2_lookup_table[t] = pow(2, t);
	}

	for (int i = 0; i < 3; i++) {
		decompressExpand_L1_Reduced_Compressed_Frame_Sparse(
			rc_fp, nx, ny, header->bit_depth,
			0, 0, 0, n_bytes_in_binary_image,
			compressedBinaryImage, deCompressedBinaryImage, compressedPixvals, deCompressedPixvals,
			pow2_lookup_table,
			frameBuffer,
			1
		);
	}

}

/*
These function are provided for testing purposes only. No external pythonic calls to these function are available. 
Decompression to only sparse format is supported for external calls.
*/
void decompressExpand_L1_Reduced_Only(FILE* fp, const char *out_fname, RCHeader *header) {}
void decompressExpand_L1_Reduced_Only_Sparse(FILE* fp, const char *out_fname, RCHeader *header) {}