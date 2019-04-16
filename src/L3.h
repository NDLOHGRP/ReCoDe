


// Used in L1
float binarize (uint8_t	 process_id,
				uint16_t *frameBuffer, 
				uint16_t *darkBuffer, 
				uint16_t epsilon, 
				DataSize h, 
				uint8_t  *binaryImage) {
	
	//uint32_t n_fg_pixels = 0;
	clock_t p_start = clock();
	
	uint16_t row, col;
	uint32_t linear_index;
	uint16_t tmp_f, tmp_d, thresh_xy;
	
	for (row = 0; row < h.ny; row++) {
		for (col = 0; col < h.nx; col++) {
			
			linear_index = row * h.nx + col;
			
			tmp_f = frameBuffer[linear_index];
			tmp_d = darkBuffer[linear_index];
			thresh_xy = tmp_d + epsilon;

			if (tmp_f > thresh_xy) {
				SetBit(binaryImage, linear_index);
				//n_fg_pixels++;
			} else {
				ClearBit(binaryImage, linear_index);
			}
		}
	}
	
	//printf("No. of foreground pixels = %lu\n", n_fg_pixels);

	clock_t p_end = clock();
	float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	return process_time;
}




float reduceFrame_L3 (  uint8_t	 process_id,
						uint16_t *frameBuffer, 
						uint16_t *darkFrame, 
						uint32_t z, 
						DataSize h, 
						uint8_t  epsilon_s,
						uint8_t  *binaryImage) {
	
	
	/*
	========================================================================================== 
	Reduce Frame to Level 1
	========================================================================================== 
	*/
	
	uint32_t frame_start_index	= z * h.nx * h.ny;
		
	float thresh_time 			= binarize (process_id, frameBuffer + frame_start_index, darkFrame, epsilon_s, h, 
											binaryImage);
	
	return thresh_time;
}



float compressFrame_L3 (uint8_t	 process_id,
						uint8_t  *binaryImage,
						DataSize h, 
						uint8_t  compression_scheme,
						uint8_t  compression_level,
						uint8_t  *compressedBinaryImage,
						uint32_t *n_compressed_bytes ) {
	
	/*
	========================================================================================== 
	Compress
	========================================================================================== 
	*/
	uint32_t n_bytes_in_binary_image = ceil((h.nx * h.ny) / 8.0);
	
	float 	 compression_time = compress_stream (compression_scheme, compression_level, 
							binaryImage, n_bytes_in_binary_image, n_compressed_bytes, compressedBinaryImage);
	
	return compression_time;
}



void reduceCompressFrame_L3 (	uint8_t	 process_id,
								uint16_t *frameBuffer, 
								uint16_t *darkFrame, 
								uint32_t frame_id,					// absolute frame index
								uint32_t z, 						// frame index relative to start of this thread
								DataSize h, 
								uint8_t  rc_operation_mode,
								uint8_t  epsilon_s,
								uint8_t  compression_scheme,
								uint8_t  compression_level,
								uint8_t  *binaryImage, 
								uint8_t  *compressedBinaryImage,
								uint32_t *n_compressed_bytes, 
								float 	 *run_metrics,
								uint8_t	 copy_compressed_frame_to_return_buffer,
								uint8_t  *compressed_frame,
								uint32_t *compressed_frame_length) {

		float reduction_time  	= reduceFrame_L3 (process_id, frameBuffer, darkFrame, z, h, 
												  epsilon_s, binaryImage);
		float compression_time  = 0.0;

		if (rc_operation_mode == RC_MODE_REDUCE_COMPRESS) {
			
			compression_time  	= compressFrame_L3 (process_id, binaryImage, 
													h, compression_scheme, compression_level,
						  							compressedBinaryImage, n_compressed_bytes);
			
		}

		uint32_t n_bytes_in_binary_image 	= ceil ((h.nx * h.ny) / 8.0);									// number of bytes needed to pack binary image
		
		/*
		========================================================================================== 
		Copy compressed data into return buffer if requested. Used in on-the-fly compression.
		========================================================================================== 
		*/
		if (copy_compressed_frame_to_return_buffer == 1) {
			
			if (rc_operation_mode 			== RC_MODE_REDUCE_ONLY) {
			
				*compressed_frame_length 	= sizeof(uint8_t) * n_bytes_in_binary_image + sizeof(uint32_t);

				memcpy(compressed_frame, &frame_id, 	sizeof(uint32_t));
				memcpy(compressed_frame, binaryImage, 	sizeof(uint8_t) * n_bytes_in_binary_image);
				
				
			} else if (rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {
				
				*compressed_frame_length 	= sizeof(uint8_t) * (*n_compressed_bytes) + 2 * sizeof(uint32_t);
				
				memcpy(compressed_frame, &frame_id, 			sizeof(uint32_t));
				memcpy(compressed_frame, n_compressed_bytes, 	sizeof(uint32_t));
				memcpy(compressed_frame, compressedBinaryImage, sizeof(uint8_t)*(*n_compressed_bytes));
				
			}
		}

		run_metrics[0] += reduction_time;
		run_metrics[1] += compression_time;
		run_metrics[2] += n_bytes_in_binary_image;
		run_metrics[3] += (*n_compressed_bytes);
		run_metrics[4] += 0;
		
		
		return;
	
}



char* reduceCompress_L3 (uint8_t	 process_id, 
						 const char  *original_filename, 
	                     const char  *out_foldername,
						 uint16_t    *frameBuffer, 				// points to the start of this thread's chunk
						 uint16_t    *darkFrame, 
						 uint32_t    frame_start_index, 
						 DataSize    h, 
						 InputParams *input_params,
						 float 		 *compression_time ) {
	
	
		
		uint32_t n_pixels_in_frame			= h.nx * h.ny;														// number of pixels in the image
		uint32_t n_bytes_in_binary_image 	= ceil (n_pixels_in_frame / 8.0);									// number of bytes needed to pack binary image


		// create data buffers
		uint8_t  *binaryImage 				= (uint8_t*) malloc (n_bytes_in_binary_image * sizeof(uint8_t));	// every bit in binaryImage has to be reset for every frame in the loop below
		uint8_t  *compressedBinaryImage 	= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
		
		// create part file
		char* part_filename = makePartFilename(process_id, original_filename, 3);
		FILE *fp = fopen (concat(out_foldername, part_filename), "wb");

	
		// serialize part header
		fwrite (&process_id, sizeof(uint8_t), 1, fp);
		
		
		uint32_t z;
		uint32_t n_compressed_bytes;

		for (z = 0; z < h.nz; z++) {

			uint32_t frame_id = frame_start_index + z;		// absolute frame index
			
			n_compressed_bytes		= 0;
			
			reduceCompressFrame_L3 (process_id, frameBuffer, darkFrame, frame_id, z, h, 
									input_params->rc_operation_mode, input_params->dark_threshold_epsilon,
									input_params->compression_scheme, input_params->compression_level, 
									binaryImage, compressedBinaryImage, &n_compressed_bytes, compression_time,
									0, NULL, NULL);
			
			if (input_params->rc_operation_mode 		== RC_MODE_REDUCE_ONLY) {

				fwrite (&frame_id, 					sizeof(uint32_t), 1, fp);
				fwrite (binaryImage , 				sizeof(uint8_t),  n_bytes_in_binary_image, fp);
				
			} else if (input_params->rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {

				fwrite (&frame_id, 					sizeof(uint32_t), 1, fp);
				fwrite (&n_compressed_bytes, 		sizeof(uint32_t), 1, fp);
				fwrite (compressedBinaryImage , 	sizeof(uint8_t),  n_compressed_bytes, fp);
				
			}
		}


		// serialize the number of frames at the end
		fwrite (&h.nz, sizeof(uint32_t), 1, fp);

		fclose (fp);


		// clean-up
		free(binaryImage);
		free(compressedBinaryImage);
		
		return part_filename;

}

void decompressExpand_L3(FILE* fp, uint16_t **frameBuffer, RCHeader **header) {
	printf("Not supported yet.");
}