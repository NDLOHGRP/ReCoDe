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
	
	char part_num[4];
	sprintf(part_num, "%03d", process_id);
	char* filename = concat(concat(original_filename, x), part_num);
	//char* filename = ".rc1_part000";
	//free(part_num);
	return filename;
}


bool startsWith(const char *str, const char *pre) {
    return strncmp(pre, str, strlen(pre)) == 0;
}

bool endsWith(const char *str, const char *suffix) {
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

void lowercase(char str[]) {
	int c = 0;
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
	char *inpfile = (char*)malloc(len+1);		// Make space for the zero.
	strncpy(inpfile, pdest, len+1);		// Copy including zero. 
	return inpfile;
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

int strcmp_ignore_case (char *a, char *b) {
	lowercase(a);
	lowercase(b);
	int c = 0;
    while (a[c] != '\0' && b[c] != '\0') {
        if (a[c] != b[c]) {
			return 0;
		}
		c++;
    }
	return 1;
}


// remove_ext: removes the "extension" from a file spec.
//   mystr is the string to process.
//   dot is the extension separator.
//   sep is the path separator (0 means to ignore).
// Returns an allocated string identical to the original but
//   with the extension removed. It must be freed when you're
//   finished with it.
// If you pass in NULL or the new string can't be allocated,
//   it returns NULL.

char *get_filename_sans_extension (char* mystr) {
	
	char dot = '.';
	char sep = PATH_SEPARATOR;
	
    char *retstr, *lastdot, *lastsep;

    // Error checks and allocate string.

    if (mystr == NULL)
        return NULL;
    if ((retstr = (char*)malloc (strlen (mystr) + 1)) == NULL)
        return NULL;

    // Make a copy and find the relevant characters.

    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, dot);
    lastsep = (sep == 0) ? NULL : strrchr (retstr, sep);

    // If it has an extension separator.

    if (lastdot != NULL) {
        // and it's before the extenstion separator.

        if (lastsep != NULL) {
            if (lastsep < lastdot) {
                // then remove it.

                *lastdot = '\0';
            }
        } else {
            // Has extension separator with no path separator.

            *lastdot = '\0';
        }
    }

    // Return the modified string.

    return retstr;
}

// checks if dir ends with a '/' and adds one if it does not
char *format_directory_path (char* dir) {
    char path_separator_str[2] = "\0";
    path_separator_str[0] = PATH_SEPARATOR;
    if (!endsWith(dir, path_separator_str)) {
        dir = concat (dir, path_separator_str);
    }
    return dir;
}

double _gettimeofday() {
	return 0.0;
}