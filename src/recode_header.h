

#define VERSION_MAJOR	0
#define VERSION_MINOR	1

#define MAX_FILENAME_LENGTH			100

#define SOURCE_FILE_TYPE_BINARY		0
#define SOURCE_FILE_TYPE_MRC		1
#define SOURCE_FILE_TYPE_OTHER		2

#define SOURCE_HEADER_POSITION_PRE		0
#define SOURCE_HEADER_POSITION_POST		1


typedef struct {
	
	/*
	identifier to demarcate an RC file. 
	Also used as process_id in part_rc files.
	*/
	uint64_t	uid;
	
	/*
	Major version number
	*/
	uint8_t		version_major;
	
	/*
	Minor version number
	*/
	uint8_t		version_minor;
	
	/*
	Data reduction to be used before compressing
	1 = Level 1, 2 = Level 2, 3 = Level 3, 4 = Level 4
	*/
	uint8_t		reduction_level;
	
	/*
	The bit depth used to store pixel intensity values. Only used for L1 and L2. Ignored in L3 and L4.
	*/
	uint8_t		bit_depth;

	/*
	Frame dimensions
	*/
	uint16_t	nx;
	uint16_t	ny;
	
	/*
	Number of frames in the dataset
	*/
	uint32_t	nz;
	
	/*
	The statistic to be stored for L2
	0 = None, 1 = Max, 2 = Mean
	For L1, L3 and L4 the following parameter is ignored
	*/
	uint8_t		L2_statistics;
	
	/*
	The centroiding scheme to be used for L4
	0 = None, 1 = Weighted Centroids, 2 = Max. Pixel, 3 = Unweighted Centroids
	For L1, L2 and L3 the parameter is ignored
	*/
	uint8_t		L4_centroiding;
	
	/*
	Compression scheme (For forward compatibility only. Currently only gzip is supported.)
	0 = none, 1 = gzip, 2 = bzip, 3 = lzma
	*/
	uint8_t		compression_scheme;
	
	/*
	Compression level (For forward compatibility only. Currently only gzip level 0 is supported.)
	gzip / bzip compression level
	*/
	uint8_t		compression_level;
	
	/*
	input file type:
	0 = Binary, 1 = MRC, 2 = Others
	used to decide the position and size of the original file's header
	if source_file_type = 1 (i.e MRC file), original file's header is pre-fixed to RC header, 
	else original file's header is post-fixed to end of compressed RC data
	*/
	uint8_t		source_file_type;
	
	/*
	input file's header size in bytes 
	used only if source_file_type is 0 or 2
	if source_file_type = 1 (i.e MRC file) this input is ignored and header size is automatically set to 1024 bytes
	*/
	uint16_t	source_header_length;
	
	/*
	position of the original source file's header in the rc file
	for compatibility MRC headers are places at the start of the file
	other headers are placed at the end
	*/
	uint8_t		source_header_position;
	
	/*
	The original source file's name. Necessary to create decompressed file.
	*/
	uint8_t		source_file_name[100];
	
	/*
	The original dark noise dataset's name. Retained for sanity purposes.
	*/
	uint8_t		dark_file_name[100];
	
	/*
	For the (ij)th pixel the threshold used for dark subtraction is estimated as:
	the maximum value of the (ij)th pixel across the z frames of the dark noise dataset + dark_threshold_epsilon
	Retained for sanity purposes.
	*/
	uint16_t	dark_threshold_epsilon;
	
	/*
	A boolean indicator: where 1 indicates the dark noise thresholds are stored as a frame, and 0 indicates they are not stored
	*/
	uint8_t		has_dark_data;
	
	/*
	Number of initial frames to skip in the original dataset
	Retained for sanity purposes.
	*/
	uint32_t	frame_offset;
	
	/*
	Number of initial frames to skip in the dark noise data set
	Retained for sanity purposes.
	*/
	uint32_t	dark_frame_offset;
	
	/*
	The bit depth used to store pixel intensity values in the original source file. 
	Retained for sanity purposes as bit depth may change in L1 and L2 and is reduced to 1 in L3 and L4.
	*/
	uint8_t		source_bit_depth;
	
	/*
	Number of frames used in dark noise data set
	Retained for sanity purposes.
	*/
	uint32_t	num_dark_frames;
	
	/*
	SHA3-256M Hash value of the uncompressed, unreduced (original) data
	256-bit hash used after decompress + expand for validation
	*/
	uint8_t		checksum[32];
	
	/*
	44 bytes reserved for future use
	to be inited to 0
	*/
	uint8_t		futures[44];
	
} RCHeader;


void create_recode_header (InputParams *input_params, uint8_t id, uint8_t *source_name, uint8_t *dark_name, RCHeader **header) {
	
	
	
	if (id != -1) {
		(*header)->uid = id;
	} else {
		(*header)->uid = 158966344846346u;
	}
	
	(*header)->version_major 			= VERSION_MAJOR;
	(*header)->version_minor 			= VERSION_MINOR;
	
	(*header)->reduction_level 			= input_params->reduction_level;
	(*header)->bit_depth 				= input_params->bit_depth;
	(*header)->nx						= input_params->num_cols;
	(*header)->ny						= input_params->num_rows;
	(*header)->nz						= input_params->num_frames;
	(*header)->L2_statistics 			= input_params->L2_statistics;
	(*header)->L4_centroiding 			= input_params->L4_centroiding;
	(*header)->compression_scheme 		= input_params->compression_scheme;
	(*header)->compression_level 		= input_params->compression_level;
	(*header)->source_file_type 		= input_params->source_file_type;
	(*header)->source_header_length 	= input_params->source_header_length;
	
	if ((*header)->source_file_type == SOURCE_FILE_TYPE_MRC) {
		(*header)->source_header_position	= SOURCE_HEADER_POSITION_PRE;
	} else {
		(*header)->source_header_position	= SOURCE_HEADER_POSITION_POST;
	}
	
	strncpy ((char*)(*header)->source_file_name, (const char*)source_name, 100);
	strncpy ((char*)(*header)->dark_file_name, (const char*)dark_name, 	 100);
	
	(*header)->dark_threshold_epsilon 	= input_params->dark_threshold_epsilon;
	(*header)->has_dark_data 			= input_params->keep_dark_data;
	(*header)->frame_offset 			= input_params->frame_offset;
	(*header)->dark_frame_offset 		= input_params->dark_frame_offset;
	(*header)->num_dark_frames 			= input_params->num_dark_frames;
	(*header)->source_bit_depth 		= input_params->source_bit_depth;
	
	int i;
	for (i = 0; i < 32; i++) {
		(*header)->checksum[i] 	= 0;
	}
	
	for (i = 0; i < 32; i++) {
		(*header)->futures[i] 	= 0;
	}
}


void parse_recode_header (FILE *fp, RCHeader **header) {
	int i;
	fread (&(*header)->uid, 					sizeof(uint64_t), 1, fp);
	fread (&(*header)->version_major, 			sizeof(uint8_t),  1, fp);
	fread (&(*header)->version_minor, 			sizeof(uint8_t),  1, fp);
	fread (&(*header)->reduction_level,			sizeof(uint8_t),  1, fp);
	fread (&(*header)->bit_depth, 				sizeof(uint8_t),  1, fp);
	fread (&(*header)->nx, 						sizeof(uint16_t), 1, fp);
	fread (&(*header)->ny, 						sizeof(uint16_t), 1, fp);
	fread (&(*header)->nz, 						sizeof(uint32_t), 1, fp);
	fread (&(*header)->L2_statistics, 			sizeof(uint8_t),  1, fp);
	fread (&(*header)->L4_centroiding, 			sizeof(uint8_t),  1, fp);
	fread (&(*header)->compression_scheme, 		sizeof(uint8_t),  1, fp);
	fread (&(*header)->compression_level, 		sizeof(uint8_t),  1, fp);
	fread (&(*header)->source_file_type, 		sizeof(uint8_t),  1, fp);
	fread (&(*header)->source_header_length, 	sizeof(uint16_t), 1, fp);
	fread (&(*header)->source_header_position,	sizeof(uint8_t),  1, fp);
	for (i = 0; i < MAX_FILENAME_LENGTH; i++) {
		fread (&(*header)->source_file_name[i],	sizeof(uint8_t),  1, fp);
	}
	for (i = 0; i < MAX_FILENAME_LENGTH; i++) {
		fread (&(*header)->dark_file_name[i],	sizeof(uint8_t),  1, fp);
	}
	fread (&(*header)->dark_threshold_epsilon, 	sizeof(uint16_t), 1, fp);
	fread (&(*header)->has_dark_data, 			sizeof(uint8_t),  1, fp);
	fread (&(*header)->frame_offset, 			sizeof(uint32_t), 1, fp);
	fread (&(*header)->dark_frame_offset, 		sizeof(uint32_t), 1, fp);
	fread (&(*header)->num_dark_frames, 		sizeof(uint32_t), 1, fp);
	fread (&(*header)->source_bit_depth, 		sizeof(uint8_t),  1, fp);
	for (i = 0; i < 32; i++) {
		fread (&(*header)->checksum[i],			sizeof(uint8_t),  1, fp);
	}
	for (i = 0; i < 44; i++) {
		fread (&(*header)->futures[i],			sizeof(uint8_t),  1, fp);
	}
}


void serialize_recode_header (FILE *fp, RCHeader *header) {
	int i;
	fwrite (&header->uid, 						sizeof(uint64_t), 1, fp);
	fwrite (&header->version_major, 			sizeof(uint8_t),  1, fp);
	fwrite (&header->version_minor, 			sizeof(uint8_t),  1, fp);
	fwrite (&header->reduction_level,			sizeof(uint8_t),  1, fp);
	fwrite (&header->bit_depth, 				sizeof(uint8_t),  1, fp);
	fwrite (&header->nx, 						sizeof(uint16_t), 1, fp);
	fwrite (&header->ny, 						sizeof(uint16_t), 1, fp);
	fwrite (&header->nz, 						sizeof(uint32_t), 1, fp);
	fwrite (&header->L2_statistics, 			sizeof(uint8_t),  1, fp);
	fwrite (&header->L4_centroiding, 			sizeof(uint8_t),  1, fp);
	fwrite (&header->compression_scheme, 		sizeof(uint8_t),  1, fp);
	fwrite (&header->compression_level, 		sizeof(uint8_t),  1, fp);
	fwrite (&header->source_file_type, 			sizeof(uint8_t),  1, fp);
	fwrite (&header->source_header_length, 		sizeof(uint16_t), 1, fp);
	fwrite (&header->source_header_position,	sizeof(uint8_t),  1, fp);
	for (i = 0; i < MAX_FILENAME_LENGTH; i++) {
		fwrite (&header->source_file_name[i],	sizeof(uint8_t),  1, fp);
	}
	for (i = 0; i < MAX_FILENAME_LENGTH; i++) {
		fwrite (&header->dark_file_name[i],		sizeof(uint8_t),  1, fp);
	}
	fwrite (&header->dark_threshold_epsilon, 	sizeof(uint16_t), 1, fp);
	fwrite (&header->has_dark_data, 			sizeof(uint8_t),  1, fp);
	fwrite (&header->frame_offset, 				sizeof(uint32_t), 1, fp);
	fwrite (&header->dark_frame_offset, 		sizeof(uint32_t), 1, fp);
	fwrite (&header->num_dark_frames, 			sizeof(uint32_t), 1, fp);
	fwrite (&header->source_bit_depth, 			sizeof(uint8_t),  1, fp);
	for (i = 0; i < 32; i++) {
		fwrite (&header->checksum[i],			sizeof(uint8_t),  1, fp);
	}
	for (i = 0; i < 44; i++) {
		fwrite (&header->futures[i],			sizeof(uint8_t),  1, fp);
	}
}


void print_recode_header (RCHeader *header) {
	int i;
	printf ("\n");
	printf ("=========================================================\n");
	printf ("ReCoDE Header\n");
	printf ("=========================================================\n");
	printf ("UID: \t\t\t %" PRIu64 "\n", header->uid);
	printf ("ReCoDE Version: \t %" PRIu8 ".", header->version_major);
	printf ("%" PRIu8 "\n", header->version_minor);
	printf ("Reduction Level: \t L%" PRIu8 " \n", header->reduction_level);
	printf ("Bit-depth: \t\t %" PRIu8 " \n", header->bit_depth);
	printf ("Number of Columns: \t %" PRIu16 " \n", header->nx);
	printf ("Number of Rows: \t %" PRIu16 " \n", header->ny);
	printf ("Number of Frames: \t %" PRIu32 " \n", header->nz);
	printf ("L2 Statistics: \t\t %" PRIu8 " \n", header->L2_statistics);
	printf ("L4 Centroiding Scheme: \t %" PRIu8 " \n", header->L4_centroiding);
	printf ("Compression: \t\t %" PRIu8 " \n", header->compression_scheme);
	printf ("Compression Level: \t %" PRIu8 " \n", header->compression_level);
	printf ("Source Type: \t\t %" PRIu8 " \n", header->source_file_type);
	printf ("Source Header Length: \t %" PRIu16 " \n", header->source_header_length);
	printf ("Source Header Position:  %" PRIu8 " \n", header->source_header_position);
	printf ("Source File Name: \t %s\n", (char*)header->source_file_name);
	printf ("Dark File Name: \t %s\n", (char*)header->dark_file_name);
	printf ("Dark Threshold Eps: \t %" PRIu16 " \n", header->dark_threshold_epsilon);
	printf ("Has Dark Data: \t\t %" PRIu8 " \n", header->has_dark_data);
	printf ("Frame Offset: \t\t %" PRIu32 " \n", header->frame_offset);
	printf ("Dark Frame Offset: \t %" PRIu32 " \n", header->dark_frame_offset);
	printf ("No. of Dark Frames: \t %" PRIu32 " \n", header->num_dark_frames);
	printf ("Source Bit-depth: \t %" PRIu8 " \n", header->source_bit_depth);
	printf ("SHA3:\t\t\t ");
	for (i = 0; i < 32; i++) {
		printf ("%" PRIu8 "", header->checksum[i]);
	}
	printf ("\n");
	printf ("=========================================================\n");
	printf ("\n");
}



void get_frame_locations_table (RCHeader *header, FILE *fp, unsigned long *frame_locations) {
	
}
