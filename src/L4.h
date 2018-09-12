



float cca (	uint16_t *foregroundImage, 
			int8_t   *clusterMap, 
			uint16_t *frameBuffer, 
			uint16_t *darkBuffer, 
			uint16_t epsilon, 
			DataSize h, 
			uint32_t n_fg_pixels, 
			float 	 *x_coor, 
			float 	 *y_coor, 
			uint32_t *n_labels) {
	

		uint32_t row, col, linear_index;

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

					do {
						//printf("1\n");
						//set the current pixel as the last element of to_be_prcessed list
						//printf("n_elem = %d\n", n_elem);
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
						//printf("3\n");

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

					//printf("Cluster %d Processed. X = %f, Y = %f.\n", label, (*x_coor)[label], (*y_coor)[label]);
				}
			}
		}

		//printf("No. of Labels: %d\n", label);
		(*n_labels) = label;

		free(bd_x);
		free(bd_y);

		return 0;
		
}



void reduceCompressFrame_L4 (	uint8_t	 process_id,
								uint16_t *frameBuffer, 
								uint16_t *darkFrame, 
								uint32_t frame_id,					// absolute frame index
								uint32_t z, 						// frame index relative to start of this thread
								DataSize h, 
								uint8_t  rc_operation_mode,
								uint16_t epsilon_s,
								uint8_t  compression_scheme,
								uint8_t  compression_level,
								float 	 *x_coor,
								float 	 *y_coor,
								uint16_t *foregroundImage,
								uint8_t  *centroidImage,
								int8_t   *foregroundTernaryMap,
								uint8_t  *compressedCentroidImage,
								uint32_t *n_compressed_bytes, 
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
		
	thresh_time = get_foreground_image (frameBuffer + frame_start_index, darkFrame, epsilon_s, h, &n_fg_pixels, foregroundImage, foregroundTernaryMap);

	
	// Connected Components Analysis to get (binary) Centroid Image
	uint32_t n_labels;
	cca_time = cca (foregroundImage, foregroundTernaryMap, frameBuffer + frame_start_index, darkFrame, 
					epsilon_s, h, n_fg_pixels, x_coor, y_coor, &n_labels);


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

	clock_t p_end = clock();
	float reduction_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;

	
	/*
	========================================================================================== 
	Compress
	========================================================================================== 
	*/
	if (rc_operation_mode == RC_MODE_REDUCE_COMPRESS) {
		compression_time = compress_stream (compression_scheme, compression_level, 
										centroidImage, n_bytes_in_binary_image, n_compressed_bytes, compressedCentroidImage);
		
	}
	
	/*
	========================================================================================== 
	Copy compressed data into return buffer if requested. Used in on-the-fly compression.
	========================================================================================== 
	*/
	if (copy_compressed_frame_to_return_buffer == 1) {
		
		if (rc_operation_mode 			== RC_MODE_REDUCE_ONLY) {
		
			*compressed_frame_length  = n_bytes_in_binary_image * sizeof(uint8_t) + sizeof(uint32_t);

			memcpy(compressed_frame,  &frame_id, 			   sizeof(uint32_t));
			memcpy(compressed_frame,  centroidImage, 		   sizeof(uint8_t) * n_bytes_in_binary_image);
			
		}  else if (rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {
			
			*compressed_frame_length  = (*n_compressed_bytes) * sizeof(uint8_t) + 2 * sizeof(uint32_t);

			memcpy(compressed_frame,  &frame_id, 			   sizeof(uint32_t));
			memcpy(compressed_frame,  n_compressed_bytes, 	   sizeof(uint32_t));
			memcpy(compressed_frame,  compressedCentroidImage, sizeof(uint8_t) * (*n_compressed_bytes));
			
		}
	}
		
	/*
	========================================================================================== 
	Log Results - to do
	========================================================================================== 
	*/
	run_metrics[0] += reduction_time;
	run_metrics[1] += compression_time;
	run_metrics[2] += n_bytes_in_binary_image;
	run_metrics[3] += (*n_compressed_bytes);
	run_metrics[4] += 0;


	return;
}


char* reduceCompress_L4 (uint8_t	 process_id, 
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
		uint16_t *foregroundImage 	 		= (uint16_t *)calloc(n_pixels_in_frame, sizeof(uint16_t));
		uint8_t  *centroidImage 			= (uint8_t *)calloc (n_bytes_in_binary_image, 1);			// allocated max space - where all pixels are fg pixels, must calloc to init to 0
		uint8_t  *foregroundTernaryMap 		= (int8_t *)calloc  (n_pixels_in_frame, sizeof(int8_t));
		uint8_t  *compressedCentroidImage	= (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
	
		// create part file
		char* part_filename = makePartFilename(process_id, original_filename, 4);
		FILE *fp = fopen (part_filename, "wb");
		
		// serialize part header
		fwrite (&process_id, sizeof(uint8_t), 1, fp);
		
		
		uint32_t z;
		uint32_t n_compressed_bytes;

		for (z = 0; z < h.nz; z++) {

			uint32_t frame_id = frame_start_index + z;		// absolute frame index
			
			n_compressed_bytes		= 0;
			
			reduceCompressFrame_L4 (process_id, frameBuffer, darkFrame, frame_id, z, h, 
									input_params->rc_operation_mode, input_params->dark_threshold_epsilon,
									input_params->compression_scheme, input_params->compression_level, 
									x_coor, y_coor, foregroundImage, centroidImage, foregroundTernaryMap, compressedCentroidImage,
									&n_compressed_bytes, compression_time,
									0, NULL, NULL);
			
			if (input_params->rc_operation_mode 		== RC_MODE_REDUCE_ONLY) {

				fwrite (&frame_id, 					sizeof(uint32_t), 1, fp);
				fwrite (centroidImage , 			sizeof(uint8_t),  n_bytes_in_binary_image, fp);
				
			} else if (input_params->rc_operation_mode 	== RC_MODE_REDUCE_COMPRESS) {

				fwrite (&frame_id, 					sizeof(uint32_t), 1, fp);
				fwrite (&n_compressed_bytes, 		sizeof(uint32_t), 1, fp);
				fwrite (compressedCentroidImage , 	sizeof(uint8_t),  n_compressed_bytes, fp);
				
			}
		}


		// serialize the number of frames at the end
		fwrite (&h.nz, sizeof(uint32_t), 1, fp);

		fclose (fp);


		// clean-up
		free(x_coor);
		free(y_coor);
		free(foregroundImage);
		free(centroidImage);
		free(foregroundTernaryMap);
		free(compressedCentroidImage);
		
		return part_filename;

}



