


float cca_with_summary_stats (	uint16_t *foregroundImage, 
								int8_t   *clusterMap, 
								uint16_t *frameBuffer, 
								uint16_t *darkBuffer, 
								uint16_t epsilon, 
								DataSize h, 
								uint32_t n_fg_pixels, 
								float 	 *x_coor, 
								float 	 *y_coor, 
								uint16_t *maxs,
								uint32_t *n_labels) {
	

		uint32_t row, col, linear_index;

        printf("Epsilon: %d\n", epsilon);
        
		//set values of cluster map boundary to zero
		/*
		for (col = 0; col < h.nx; col++) {
			clusterMap[col * h.ny] = 0;
			clusterMap[col * h.ny + (h.ny - 1)] = 0;
		}

		for (row = 0; row < h.ny; row++) {
			clusterMap[row * h.nx] = 0;
			clusterMap[row * h.nx + (h.nx - 1)] = 0;
		}
		*/

		uint32_t x, y;

		
		uint32_t l0, l1, y0, y1, x0, x1;
		for (x = 0; x < h.nx; x++) {
			y0 = 0;
			y1 = h.ny - 1;
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
		

		uint32_t *bd_x = (uint32_t *) malloc(n_fg_pixels * sizeof(uint32_t));
		uint32_t *bd_y = (uint32_t *) malloc(n_fg_pixels * sizeof(uint32_t));

		
		uint32_t label = 0;
		uint32_t temp = 0;


		int nb[8][2] = { { -1,0 }, { -1,-1 }, { 0,-1 }, { 1,-1 }, { 1,0 }, { 1,1 }, { 0,1 }, { -1,1 } };
		int MAX_NEIGHBORS = 8;

		//Identify and label clusters
		for (y = 0; y < h.ny; y++) {
		
			for (x = 0; x < h.nx; x++) {

				if (clusterMap [y * h.nx + x] == -1) {				// found a fg pixel - start labelling
				
					int32_t n_elem = 0;								// num elements left in the to_be_processed_queue
					uint32_t intensity_sum = 0; 

					uint32_t n, x_curr, y_curr, x_neighbor, y_neighbor;
				
					
					bd_x[n_elem] = x;								// set current pixels to: to_be_processed list
					bd_y[n_elem] = y;
					
					label++;										// create new label 
					//printf("New Cluster Label: %d\n", label);

					// malloced array, so initialize
					x_coor[label] = 0;
					y_coor[label] = 0;
                    maxs[label] = (uint16_t)(floor(pow(2.0, (double)h.dtype))-1);
                    //printf("Max Possible = %" PRIu16 "\n", maxs[label]);

					do {
						
						x_curr = bd_x[n_elem];
						y_curr = bd_y[n_elem];

						// assign current pixel the current lable
						linear_index = y_curr * h.nx + x_curr;
						clusterMap [linear_index] = 2;					// label it, so that it is not revisited
						
						// sum of pixel intensities (denominator for centroid computation)
						intensity_sum += foregroundImage[linear_index];

						// pixel intensity weighted positions (numerator for centroid computation)
						x_coor[label] += x_curr * foregroundImage[linear_index];
						y_coor[label] += y_curr * foregroundImage[linear_index];
						if (maxs[label] > foregroundImage[linear_index]) {
							maxs[label] = foregroundImage[linear_index];
						}

						// decrease queue counter
						n_elem--;


						// check current pixel's neighbors
						for (n = 0; n < MAX_NEIGHBORS; n++)
						{
							x_neighbor = x_curr + nb[n][0];
							y_neighbor = y_curr + nb[n][1];
							linear_index = y_neighbor * h.nx + x_neighbor;
							

							if (clusterMap [linear_index] == -1)		// if neighbor is fg pixel, add it to the queue
							{
								n_elem++;
								bd_x[n_elem] = x_neighbor;
								bd_y[n_elem] = y_neighbor;
								clusterMap [linear_index] = 2;			// also label them, so that they are not revisited
								
							}

						}

					} while (n_elem >= 0);		// queue is empty => this cluster is done


					// compute the centroids
					x_coor[label] /= (intensity_sum * 1.0);
					y_coor[label] /= (intensity_sum * 1.0);

				}
			}
		}

		printf("No. of Labels: %d\n", label);
		(*n_labels) = label;

		free(bd_x);
		free(bd_y);

		return 0;
		
}


float pack_summary_stats (	uint8_t	 process_id,
							uint16_t *eventSummaryStats,
							uint32_t n_values,
							uint8_t  bit_depth, 
							uint8_t  *packedSummaryStats) {
	
	int p, n, linear_index, nth_bit_of_pth;

	for (p=0; p<n_values; p++) {
		
		// Assumes LITTLE-ENDIAN Byte Order
		for (n=0; n<bit_depth; n++) {
			nth_bit_of_pth = (eventSummaryStats[p] & ( 1 << n )) >> n;
			if (nth_bit_of_pth == 1) {
				linear_index = p*bit_depth + n;
				SetBit(packedSummaryStats, linear_index);
			}
		}
	}
	
	//printf("MAX_VAL for given bit_depth: %hu\n", MAX_VAL);
	
}


void reduceCompressFrame_L2 (	uint8_t	 process_id,
								uint16_t *frameBuffer, 
								uint16_t *darkFrame, 
								uint32_t frame_id,					// absolute frame index
								uint32_t z, 						// frame index relative to start of this thread
								DataSize h, 
								uint8_t  rc_operation_mode,
								uint16_t epsilon_s,
								uint8_t	 bit_depth,
								uint8_t  compression_scheme,
								uint8_t  compression_level,
								float 	 *x_coor,
								float 	 *y_coor,
								uint16_t *foregroundImage,
								uint8_t  *centroidImage,
								uint16_t *eventSummaryStats,
								int8_t   *foregroundTernaryMap,
								uint8_t  *compressedCentroidImage,
								uint8_t  *packedSummaryStats,
								uint8_t  *compressedSummaryStats,
								uint32_t *n_bytes_in_summary_stats, 
								uint32_t *n_compressed_bytes_centroid_image, 
								uint32_t *n_compressed_bytes_summary_stats, 
								float 	 *run_metrics,
								uint8_t	 copy_compressed_frame_to_return_buffer,
								uint8_t  *compressed_frame,
								uint32_t *compressed_frame_length) {

							
	float thresh_time, cca_time;
	float compression_time = 0;

	clock_t p_start = clock();
	
	uint32_t n_bytes_in_binary_image = ceil((h.nx * h.ny) / 8.0);

	/*
	========================================================================================== 
	Reduce
	========================================================================================== 
	*/
	
	// Dark Subtraction
	uint32_t n_fg_pixels = 0;
	uint32_t frame_start_index = frame_id * h.nx * h.ny;
		
	thresh_time = get_foreground_image (frameBuffer + frame_start_index, darkFrame, epsilon_s, h, 
										&n_fg_pixels, foregroundImage, foregroundTernaryMap);

    printf ("No. Fg. Pixels = %" PRIu32 "\n", n_fg_pixels);
	
	// Connected Components Analysis to get (binary) Centroid Image
	uint32_t n_labels;
	cca_time = cca_with_summary_stats ( foregroundImage, foregroundTernaryMap, frameBuffer + frame_start_index, darkFrame, 
										epsilon_s, h, n_fg_pixels, x_coor, y_coor, eventSummaryStats, &n_labels);

    /*
    int p;
    for (p=0; p<n_labels; p++) {
        printf("X = %f, Y = %f, Max = %" PRIu16 "\n", x_coor[p], y_coor[p], eventSummaryStats[p]);
    }
    serializeFrames (foregroundImage, h.nx, h.ny, "foregroundImage_0.txt", 0);
    */
                                        
	// Make Binary / Centroid Image
	uint32_t i, linear_index;
	uint32_t n_pixels_in_frame = h.nx * h.ny;
	for (i = 0; i < n_pixels_in_frame; i++) {
		ClearBit(centroidImage, i);
	}
	for (i = 1; i <= n_labels; i++) {
		linear_index = floor(y_coor[i]) * h.nx + floor(x_coor[i]);
		SetBit(centroidImage, linear_index);
	}

	// pack summary stats
	*n_bytes_in_summary_stats = n_labels * ceil(bit_depth/8.0);
	pack_summary_stats (process_id, eventSummaryStats, n_labels, bit_depth, packedSummaryStats);


	clock_t p_end = clock();
	float reduction_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;
	
	/*
	========================================================================================== 
	Compress
	========================================================================================== 
	*/
	if (rc_operation_mode == RC_MODE_REDUCE_COMPRESS) {
		compression_time = compress_stream (compression_scheme, compression_level, 
					centroidImage, n_bytes_in_binary_image, n_compressed_bytes_centroid_image, compressedCentroidImage);
													
		compression_time += compress_stream (compression_scheme, compression_level, 
					packedSummaryStats, *n_bytes_in_summary_stats, n_compressed_bytes_summary_stats, compressedSummaryStats);
		
	}
	
	/*
	========================================================================================== 
	Copy compressed data into return buffer if requested. Used in on-the-fly compression.
	========================================================================================== 
	*/
	if (copy_compressed_frame_to_return_buffer == 1) {
		
		if (rc_operation_mode 			== RC_MODE_REDUCE_ONLY) {
		
			*compressed_frame_length  = n_bytes_in_binary_image * sizeof(uint8_t) + 
										(*n_bytes_in_summary_stats) * sizeof(uint8_t) + 
										sizeof(uint32_t);

			memcpy(compressed_frame,  &frame_id, 		  sizeof(uint32_t));
			memcpy(compressed_frame,  centroidImage, 	  sizeof(uint8_t) * n_bytes_in_binary_image);
			memcpy(compressed_frame,  packedSummaryStats, sizeof(uint8_t) * (*n_bytes_in_summary_stats));
			
		}  else if (rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {
			
			*compressed_frame_length  = (*n_compressed_bytes_centroid_image) * sizeof(uint8_t) + 
										(*n_compressed_bytes_summary_stats) * sizeof(uint8_t) + 
										3 * sizeof(uint32_t);

			memcpy(compressed_frame,  &frame_id, 			   			 sizeof(uint32_t));
			memcpy(compressed_frame,  n_compressed_bytes_centroid_image, sizeof(uint32_t));
			memcpy(compressed_frame,  n_compressed_bytes_summary_stats,  sizeof(uint32_t));
			memcpy(compressed_frame,  compressedCentroidImage, 			 sizeof(uint8_t) * (*n_compressed_bytes_centroid_image));
			memcpy(compressed_frame,  packedSummaryStats, 				 sizeof(uint8_t) * (*n_compressed_bytes_summary_stats));
			
		}
	}
		
	/*
	========================================================================================== 
	Log Results - to do
	========================================================================================== 
	*/
	run_metrics[0] += reduction_time;
	run_metrics[1] += compression_time;
	run_metrics[2] += (n_bytes_in_binary_image + *n_bytes_in_summary_stats);
	run_metrics[3] += (*n_compressed_bytes_centroid_image + *n_compressed_bytes_summary_stats);
	run_metrics[4] += 0;


	return;
}


char* reduceCompress_L2 (uint8_t	 process_id, 
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
		float 	 *x_coor 			 		= (float *) malloc  (n_pixels_in_frame * sizeof(float));
		float 	 *y_coor 			 		= (float *) malloc  (n_pixels_in_frame * sizeof(float));
		uint16_t *foregroundImage 	 		= (uint16_t *)calloc(n_pixels_in_frame,  sizeof(uint16_t));
		uint8_t  *centroidImage 			= (uint8_t *)calloc (n_bytes_in_binary_image, 1);
		uint16_t *eventSummaryStats			= (uint16_t *)calloc(n_pixels_in_frame,  sizeof(uint16_t));
		uint8_t  *foregroundTernaryMap 		= (int8_t *)calloc  (n_pixels_in_frame,  sizeof(int8_t));
		uint8_t  *compressedCentroidImage	= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
		uint8_t  *packedSummaryStats		= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
		uint8_t  *compressedSummaryStats	= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
	
		// create part file
		char* part_filename = makePartFilename(process_id, original_filename, 2);
		FILE *fp = fopen (part_filename, "wb");
		
		// serialize part header
		fwrite (&process_id, sizeof(uint8_t), 1, fp);
		
		
		uint32_t z;
		uint32_t n_bytes_in_summary_stats;
		uint32_t n_compressed_bytes_centroid_image;
		uint32_t n_compressed_bytes_summary_stats;

		for (z = 0; z < h.nz; z++) {

			uint32_t frame_id = frame_start_index + z;		// absolute frame index
			
			n_bytes_in_summary_stats 			= 0;
			n_compressed_bytes_centroid_image	= 0;
			n_compressed_bytes_summary_stats	= 0;
			
			reduceCompressFrame_L2 (process_id, frameBuffer, darkFrame, frame_id, z, h, 
									input_params->rc_operation_mode, input_params->dark_threshold_epsilon, input_params->bit_depth,
									input_params->compression_scheme, input_params->compression_level, 
									x_coor, y_coor, foregroundImage, centroidImage, eventSummaryStats,
									foregroundTernaryMap, compressedCentroidImage, packedSummaryStats, compressedSummaryStats,
									&n_bytes_in_summary_stats, &n_compressed_bytes_centroid_image, &n_compressed_bytes_summary_stats, 
									compression_time,
									0, NULL, NULL);
			
			if (input_params->rc_operation_mode 		== RC_MODE_REDUCE_ONLY) {

				fwrite (&frame_id, 			sizeof(uint32_t), 1, fp);
				fwrite (centroidImage, 		sizeof(uint8_t),  n_bytes_in_binary_image, fp);
				fwrite (eventSummaryStats, 	sizeof(uint8_t),  n_bytes_in_summary_stats, fp);
				
			} else if (input_params->rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {

				fwrite (&frame_id, 							sizeof(uint32_t), 1, fp);
				fwrite (&n_compressed_bytes_centroid_image, sizeof(uint32_t), 1, fp);
				fwrite (&n_compressed_bytes_summary_stats, 	sizeof(uint32_t), 1, fp);
				fwrite (compressedCentroidImage, 			sizeof(uint8_t),  n_compressed_bytes_centroid_image, fp);
				fwrite (compressedSummaryStats, 			sizeof(uint8_t),  n_compressed_bytes_summary_stats, fp);
				
			}
		}


		// serialize the number of frames at the end
		fwrite (&h.nz, sizeof(uint32_t), 1, fp);

		fclose (fp);

		printf("Here...");

		// clean-up
		free(x_coor);						printf("1, ");
		free(y_coor);						printf("2, ");
		free(foregroundImage);				printf("3, ");
		free(centroidImage);				printf("4, ");
		free(eventSummaryStats);			printf("5, ");
		free(foregroundTernaryMap);			printf("6, ");
		free(compressedCentroidImage);		printf("7, ");
		free(packedSummaryStats);			printf("8, ");
		free(compressedSummaryStats);		printf("9..");
		
		printf("but not here...\n");
		
		return part_filename;

}



