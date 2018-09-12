char* concat(const char *s1, const char *s2) {
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = (char*) malloc(len1+len2+1);//+1 for the null-terminator
    //check for errors in malloc
    memcpy(result, s1, len1);
    memcpy(result+len1, s2, len2+1);//+1 to copy the null-terminator
    return result;
}


char* makePartFilename(unsigned char process_id, const char *original_filename, unsigned int reduction_level) {
	
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
	
	char* part_num = malloc(sizeof(char)*3);
	sprintf(part_num, "%03d", process_id);
	char* filename = concat(concat(original_filename, x), part_num);
	free(part_num);
	return filename;
}


bool startsWith(const char *str, const char *pre) {
    return strncmp(pre, str, strlen(pre)) == 0;
}

void lowercase(char str[]) {
	int c;
	while (str[c] != '\0') {
		if (str[c] >= 'A' && str[c] <= 'Z') {
			str[c] = str[c] + 32;
		}
		c++;
	}
}

char* getFilenameFromPath (char *filepath) {

	// Search backwards for last backslash in filepath 
	char *pdest = strrchr(filepath, PATH_SEPARATOR);

	// if backslash not found in filepath
	if(pdest == NULL ) {
		pdest = filepath; 	// The whole name is a file in current path? 
	} else {
		pdest++; 			// Skip the backslash itself.
	}

	// extract filename from file path
	int len = strlen(pdest);
	char *inpfile = malloc(len+1);		// Make space for the zero.
	strncpy(inpfile, pdest, len+1);		// Copy including zero. 
	return inpfile;
}

const char* get_filename_sans_extension (const char* filename) {
	
}

char *trim(char *str) {
	
	char *end;

	// Trim leading space
	while(isspace((unsigned char)*str)) str++;

	if(*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace((unsigned char)*end)) end--;

	// Write new null terminator
	*(end+1) = 0;

	return str;
}