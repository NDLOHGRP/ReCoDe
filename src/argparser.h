
typedef struct {
	
	unsigned char	reduction_level;
	
	/*
	rc_operation_mode can be 0 => reduction only or 1 => reduction and compression, 
	there is also a compression only mode, rc_operation_mode = 2 but this is used only for diagnostics
	Default is rc_operation_mode = 1.
	*/
	unsigned char	rc_operation_mode;
	
	/*
	For the (ij)th pixel the threshold used for dark subtraction is estimated as:
	the maximum value of the (ij)th pixel across the z frames of the dark noise dataset + dark_threshold_epsilon
	*/
	unsigned short	dark_threshold_epsilon;
	
	unsigned char	bit_depth;
	
	/*
	The bit depth of pixel intensity values in the source file.
	If source file is MRC the bit depth as per MRC header is used
	If source file is not MRC this value must be explicitly provided
	*/
	unsigned char	source_bit_depth;
	
	/*
	Frame dimensions
	If source file is MRC, these parameters are read from MRC header, else user provided values are used
	*/
	unsigned short	num_cols;
	unsigned short	num_rows;
	
	/*
	Number of frames to extract from dataset i
	If source file is MRC and the user provided value is -1, the number of frames as per MRC header is used (i.e. all frames are used)
	If source file is MRC and the user provided value is not -1, the user specified number of frames are read
	If source file is not MRC this value has to be explicitly provided
	*/
	unsigned long	num_frames;
	
	/*
	Number of initial frames to skip in the original dataset
	*/
	unsigned long	frame_offset;
	
	/*
	Number of frames to be used for dark noise estimation
	If source file is MRC and the user provided value is -1, the number of frames as per MRC header is used (i.e. all frames are used)
	If source file is MRC and the user provided value is not -1, the user specified number of frames are read
	If source file is not MRC this value has to be explicitly provided
	*/
	unsigned long	num_dark_frames;
	
	/*
	Number of initial frames to skip in the dark noise data set
	*/
	unsigned long	dark_frame_offset;
	
	unsigned char	keep_part_files;
	unsigned char	num_threads;
	unsigned char	L2_statistics;		// 0 = None (default for L1, L3 and L4), 1 = Max, 2 = Mean
	unsigned char	L4_centroiding;		// 0 = None (default for L1, L2 and L3), 1 = Weighted Centroids, 1 = Max. Pixel, 2 = Unweighted Centroids
	unsigned char	compression_scheme;	// 0 = none, 1 = gzip, 2 = bzip, 3 = lzma (For forward compatibility only. Currently only gzip is supported. The parameter value will be ignored and automatically set to 1.)
	unsigned char	compression_level;	// gzip / bzip compression level (For forward compatibility only. Currently only gzip level 0 is supported. The parameter value will be ignored and automatically set to 0.)
	unsigned char	source_file_type;
	unsigned short	source_header_length;
	
	unsigned char	keep_dark_data;
	
} InputParams;



get_input_params (const char* argfilename, InputParams **params) {
	
	FILE *fp = fopen ( argfilename, "r" );
	if ( fp != NULL ) {
		
		char line [ 1024 ];
		char *name;
		int value;
		char *search = "=";

		while ( fgets ( line, sizeof line, fp ) != NULL ) {
		
			if (!startsWith(line, "#") && strcmp(trim(line), "") != 0) {
				
				name = trim(strtok(line, search));
				lowercase(name);
				value = atoi(trim(strtok(NULL, search)));
				
				if (strcmp(name, "reduction_level") == 0) {
					
					(*params)->reduction_level = (unsigned char)value;
					
				} else if (strcmp(name, "rc_operation_mode") == 0) {
					
					(*params)->rc_operation_mode = (unsigned short)value;
					
				} else if (strcmp(name, "dark_threshold_epsilon") == 0) {
					
					(*params)->dark_threshold_epsilon = (unsigned short)value;
					
				} else if (strcmp(name, "bit_depth") == 0) {
					
					(*params)->bit_depth = (unsigned char)value;
					
				} else if (strcmp(name, "source_bit_depth") == 0) {
					
					(*params)->source_bit_depth = (unsigned char)value;
					
				} else if (strcmp(name, "num_cols") == 0) {
					
					(*params)->num_cols = (unsigned short)value;
					
				} else if (strcmp(name, "num_rows") == 0) {
					
					(*params)->num_rows = (unsigned short)value;
					
				} else if (strcmp(name, "num_frames") == 0) {
					
					(*params)->num_frames = (signed long)value;
					
				} else if (strcmp(name, "frame_offset") == 0) {
					
					(*params)->frame_offset = (signed long)value;
					
				} else if (strcmp(name, "num_dark_frames") == 0) {
					
					(*params)->num_dark_frames = (signed long)value;
					
				} else if (strcmp(name, "dark_frame_offset") == 0) {
					
					(*params)->dark_frame_offset = (signed long)value;
					
				} else if (strcmp(name, "keep_part_files") == 0) {
					
					(*params)->keep_part_files = (unsigned char)value;
					
				} else if (strcmp(name, "num_threads") == 0) {
					
					(*params)->num_threads = (unsigned char)value;
					
				} else if (strcmp(name, "l2_statistics") == 0) {
					
					(*params)->L2_statistics = (unsigned char)value;
					
				} else if (strcmp(name, "l4_centroiding") == 0) {
					
					(*params)->L4_centroiding = (unsigned char)value;
					
				} else if (strcmp(name, "compression_scheme") == 0) {
					
					(*params)->compression_scheme = (unsigned char)value;
					
				} else if (strcmp(name, "compression_level") == 0) {
					
					(*params)->compression_level = (unsigned char)value;
					
				} else if (strcmp(name, "source_file_type") == 0) {
					
					(*params)->source_file_type = (unsigned char)value;
					
				} else if (strcmp(name, "source_header_length") == 0) {
					
					(*params)->source_header_length = (unsigned short)value;
					
				} else if (strcmp(name, "keep_dark_data") == 0) {
					
					(*params)->keep_dark_data = (unsigned char)value;
					
				} else {
					
					printf("Unknown Argument: %s\n", name);
					
				}

				
			}
		}
		
		fclose ( fp );
		
	} else {
		
		(*params)->reduction_level 			= 1;
		(*params)->rc_operation_mode 		= 1;
		(*params)->dark_threshold_epsilon	= 0;
		(*params)->bit_depth 				= 12;
		(*params)->source_bit_depth 		= 12;
		(*params)->num_cols 				= 4096;
		(*params)->num_rows 				= 4096;
		(*params)->num_frames 				= -1;
		(*params)->frame_offset				= 0;
		(*params)->num_dark_frames 			= -1;
		(*params)->dark_frame_offset		= 0;
		(*params)->keep_part_files 			= 0;
		(*params)->num_threads 				= 1;
		(*params)->L2_statistics 			= 0;
		(*params)->L4_centroiding 			= 0;
		(*params)->compression_scheme 		= 0;
		(*params)->compression_level 		= 0;
		(*params)->source_file_type 		= 0;
		(*params)->source_header_length 	= 0;
	}
	
}

compile_missing_params (MRCHeader *mrc_header, InputParams **params) {
	
	if ((*params)->source_file_type == 1) {
		(*params)->num_cols 				= mrc_header->nx;
		(*params)->num_rows 				= mrc_header->ny;
		(*params)->source_bit_depth 		= mrc_header->bit_depth;
		if ((*params)->num_frames == -1) {
			(*params)->num_frames			= mrc_header->nz;
		}
		(*params)->source_header_length 	= 1024;
	}
	
}


print_params (InputParams *params) {
	printf("Reduction Level: %d\n", params->reduction_level);
	printf("Num. Dark Frames: %d\n", params->num_dark_frames);
	printf("Dark Frame Offset: %d\n", params->dark_frame_offset);
	printf("Num. Frames: %d\n", params->num_frames);
	printf("Frame Offset: %d\n", params->frame_offset);
}


int validate_input_params (InputParams *input_params, uint8_t *source_name, uint8_t *dark_name) {
	
	if (strlen(source_name) > (MAX_FILENAME_LENGTH-1)) {
		printf("Source filename must be less than 100 characters long.\n");
		printf("Unable to continue.\n");
		exit(0);
	}
	
	if (strlen(dark_name) > (MAX_FILENAME_LENGTH-1)) {
		printf("Dark filename must be less than 100 characters long.\n");
		printf("Unable to continue.\n");
		exit(0);
	}
	
}
