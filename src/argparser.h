
void get_input_params (const char* argfilename, InputParams **params) {
	
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
					
				} else if (strcmp(name, "target_bit_depth") == 0) {
					
					(*params)->target_bit_depth = (unsigned char)value;
					
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
					
				} else if (strcmp(name, "keep_dark_data") == 0) {
					
					(*params)->keep_dark_data = (unsigned char)value;
					
				} else if (strcmp(name, "source_file_type") == 0) {

					(*params)->source_file_type = (unsigned char)value;

				} else if (strcmp(name, "source_header_length") == 0) {

					(*params)->source_header_length = (unsigned short)value;

				} else if (strcmp(name, "dark_file_type") == 0) {

					(*params)->dark_file_type = (unsigned char)value;

				} else if (strcmp(name, "dark_header_length") == 0) {

					(*params)->dark_header_length = (unsigned short)value;

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
		(*params)->target_bit_depth 		= 12;
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
		(*params)->keep_dark_data			= 0;
		(*params)->source_file_type 		= 0;
		(*params)->source_header_length 	= 0;
		(*params)->dark_file_type			= 0;
		(*params)->dark_header_length		= 0;
	}
	
}

void compile_missing_params (MRCHeader *mrc_header, InputParams **params) {
	
	if ((*params)->source_file_type == SOURCE_FILE_TYPE_MRC) {
		(*params)->num_cols 				= mrc_header->dimensions[0];
		(*params)->num_rows 				= mrc_header->dimensions[1];
		(*params)->source_bit_depth 		= getBitDepth (mrc_header);
		if ((*params)->num_frames == -1) {
			(*params)->num_frames			= mrc_header->dimensions[2];
		}
		(*params)->source_header_length 	= 1024;
	}
	
	if ((*params)->dark_file_type == SOURCE_FILE_TYPE_MRC) {
		(*params)->dark_header_length = 1024;
	}
}


void print_params (InputParams *params) {
	printf("Reduction Level: %d\n", params->reduction_level);
	printf("Num. Dark Frames: %d\n", params->num_dark_frames);
	printf("Dark Frame Offset: %d\n", params->dark_frame_offset);
	printf("Num. Frames: %d\n", params->num_frames);
	printf("Frame Offset: %d\n", params->frame_offset);
}


int validate_input_params (InputParams *input_params, uint8_t *source_name, uint8_t *dark_name) {
	
	if (strlen((const char*)source_name) > (MAX_FILENAME_LENGTH-1)) {
		printf("Source filename must be less than 100 characters long.\n");
		printf("Unable to continue.\n");
		exit(0);
	}
	
	if (strlen((const char*)dark_name) > (MAX_FILENAME_LENGTH-1)) {
		printf("Dark filename must be less than 100 characters long.\n");
		printf("Unable to continue.\n");
		exit(0);
	}
	
}
