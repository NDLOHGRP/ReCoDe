int _sequence_file_read_test () {

	/*====================Testing========================*/
	//const char* filename = "D:/cbis/images/20161207/14-37-50.811.seq";
	//const char* filename = "Z:/26-April-2019/Live_Acquisition_1/RamDisk_LeftOver/15-10-17.540.seq";
	const char* filename = "D:/cbis/images/SequenceBlocks/17-21-33.138_Dark_Ref_3.seq";
	DataSize d = { 4096, 512, 2000, 12 };
	uint64_t sz_darkBuffer = d.nx * d.ny * d.nz * sizeof(uint16_t);
	uint16_t *darkBuffer = (uint16_t *)malloc(sz_darkBuffer);
	
	SEQHeader *seq_header = (SEQHeader *)malloc(sizeof(SEQHeader));
	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		printf("Failed to open file: '%s' for reading.", filename);
		exit(1);
	}
	getSEQHeader(fp, &seq_header);
	uint64_t frame_offset = 0;
	uint64_t available_frames = getSEQFrames(0, fp, darkBuffer, frame_offset, d.nz, seq_header);
	printf("Available Frames = %" PRIu64 "\n", available_frames);
	uint64_t i, j, frame_start_index, linear_index;
	uint16_t max_val = 0;
	int x[] = { 99, 1111, 2222, 3333, 4000 };
	int y[] = { 110, 220, 330, 440, 500 };
	for (i = 0; i <  available_frames; i++) {
		frame_start_index = i*d.nx*d.ny;
		printf("Frame %" PRIu64 ": ", i + frame_offset);
		for (j = 0; j < 5; j++) {
			linear_index = y[j] * d.nx + x[j];
			printf("%d,%d=%" PRIu16 ", ", y[j], x[j], darkBuffer[frame_start_index + linear_index]);
		}
		printf("\n");
	}
	printf("Max = %" PRIu16 "\n", max_val);
	fclose(fp);

	exit(0);

	/*
	createSEQFiles(darkBuffer, seq_header, 50, "D:/cbis/GitHub/ReCoDe/scratch/temp/", 10);
	return 0;
	*/
}