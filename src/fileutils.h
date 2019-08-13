#define MAX_PART_FILE_NAME_LENGTH	100

void writeToFile (const char *filename, uint8_t *data, int dataLength) {
	
	FILE * fp;
	fp = fopen (filename, "w+");
	
	int i = 0;
	for (i = 0; i < dataLength; i++) {
		fprintf(fp, "%u ", data[i]);
	}
	
	fclose(fp);
	
}


void serializeFrames (uint16_t *data, uint32_t nx, uint32_t ny, const char* filename, int nFrames) {
	
	FILE *fp = fopen (filename, "w+");
	
	uint32_t row, col, linear_index;
	for (row = 0; row < ny; row++) {
		for (col = 0; col < nx; col++) {
			linear_index = row * nx + col;
			fprintf(fp, "%u ", data[linear_index]);
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
}


char** getPartFilenames (const char* foldername, const char* basefilename, int reduction_level, int* num_partfiles) {
	
	const char* x;
	if (reduction_level == 1) {
		x = ".rc1_part";
	} else if (reduction_level == 2) {
		x = ".rc2_part";
	} else if (reduction_level == 3) {
		x = ".rc3_part";
	}  else if (reduction_level == 4) {
		x = ".rc4_part";
	} else {
		x = ".rcX_part";
	}
	
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (foldername)) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			const char* fname = concat(basefilename, x);
			if (startsWith(ent->d_name, fname)) {
				(*num_partfiles)++;
			}
		}
		closedir (dir);
	}
	
	int count = 0;
	char** part_filenames = (char**)malloc((*num_partfiles)*sizeof(char*));
	if ((dir = opendir (foldername)) != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			const char* fname = concat(basefilename, x);
			if (startsWith(ent->d_name, fname)) {
				part_filenames[count] = (char*)malloc(MAX_PART_FILE_NAME_LENGTH*sizeof(char));
				sprintf(part_filenames[count++], ent->d_name);
			}
		}
		closedir (dir);
	}
	
	return part_filenames;
}

void saveFrame (const char* filename, unsigned short* buffer, DataSize d, int framenum) {
	
	int r = 1;
	int c = 1;
	
	FILE * fp;

	fp = fopen (filename, "w+");
	
	long int offset = d.ny*d.nx*framenum;
	for (r = 0; r < d.ny; r++) {
		for (c = 0; c < d.nx; c++) {
			fprintf(fp, "%d ", buffer[offset + r*d.nx + c]);
		}
		fprintf(fp, "\n");
	}
	
	fclose(fp);
	
}

