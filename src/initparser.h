
#pragma warning(disable:4996)
// define the ID values to indentify the option
enum { OPT_HELP, OPT_RC, OPT_DE, ARG_VERB, ARG_IMG, ARG_DRK, ARG_PRM, ARG_DIR, ARG_LOG, ARG_NAM, ARG_VAF };

CSimpleOpt::SOption g_rgOptions[] = {
	{ OPT_RC,	 _T("-rc"),    SO_NONE },	  // "-rc: reduce-compress"
	{ OPT_DE,	 _T("-de"),    SO_NONE },	  // "-de: decompress-expand"
	{ ARG_VERB,  _T("-v"),     SO_REQ_SEP },  // "-v verbosity"
	{ ARG_IMG,   _T("-i"),     SO_REQ_SEP },  // "-i imagefile; -i recode_file"
	{ ARG_DRK,   _T("-d"),     SO_REQ_SEP },  // "-d darkfile"
	{ ARG_PRM,   _T("-p"),     SO_REQ_SEP },  // "-p params_file"
	{ ARG_DIR,   _T("-o"),     SO_REQ_SEP },  // "-o output_directory"
	{ ARG_LOG,   _T("-l"),     SO_REQ_SEP },  // "-l log_file"
	{ ARG_NAM,   _T("-n"),     SO_REQ_SEP },  // "-n name for the run"
	{ ARG_VAF,   _T("-vf"),    SO_REQ_SEP },  // "-vf gap between validation frames"
	{ OPT_HELP,  _T("-?"),     SO_NONE },	  // "-?"
	{ OPT_HELP,  _T("-h"),     SO_NONE },	  // "-h"
	{ OPT_HELP,  _T("-help"),  SO_NONE },	  // "-help"
	{ OPT_HELP,  _T("--help"), SO_NONE },	  // "--help"
	SO_END_OF_OPTIONS
};

typedef struct {
	int8_t mode;
	uint8_t verbosity;
	int64_t validation_frame_gap;
	char* image_filename;
	char* dark_filename;
	char* params_filename;
	char* output_directory;
	char* log_filename;
	char* run_name;
} InitParams;

int char_to_int(char *d) {
	int val = atoi(d);
	return (int)val;
}

int wchar_to_int(wchar_t *src) {
	uint64_t L = wcslen(src);
	char *dst = (char*)malloc(sizeof(char)*L);
	wcstombs(dst, src, L);
	int val = atoi(dst);
	return (int)val;
}

char* wchar_to_char(wchar_t *src) {
	uint64_t L = wcslen(src);
	char *dst = (char*)malloc(sizeof(char)*L + 1);
	wcstombs(dst, src, L);
	dst[L] = '\0';
	return dst;
}

/*
int wchar_to_int(char *src) {
	uint64_t L = wcslen(src);
	char *dst = (char*)malloc(sizeof(char)*L);
	wcstombs(dst, src, L);
	int val = atoi(dst);
	return (int)val;
}

char* wchar_to_char(wchar_t *src) {
	uint64_t L = wcslen(src);
	char *dst = (char*)malloc(sizeof(char)*L + 1);
	wcstombs(dst, src, L);
	dst[L] = '\0';
	return dst;
}
*/

// show the usage of this program
void ShowUsage() {
	printf("Usage:\n");
	printf("ReCoDe -rc -i ARG -d ARG -p ARG -o ARG [-v ARG] [-?] [--help]\n");
	printf("ReCoDe -de -i ARG -o ARG [-v ARG] [-?] [--help]\n");
	printf("\n");
	printf("-rc:    Perform Reduction-Compression (Either -rc or -de must be specified)\n");
	printf("-de:    Perform Decompression-Expansion (Either -rc or -de must be specified)\n");
	printf("-i:     (Required) Image file to be compressed when using -rc and ReCoDe file to be decompressed when using -de\n");
	printf("-o:     (Required) Output directory\n");
	printf("-d:     Dark file (Required when using -rc)\n");
	printf("-p:     Params file (Required when using -rc)\n");
	printf("-v:     Verbosity level (0 or 1, Optional)\n");
	printf("-l:     Log file name (Optional)\n");
	printf("-n:     Run name (Optional). Used while logging.\n");
	printf("-vf:    Gap between validation frames. (Optional). If not specified no validation frames are saved.\n");
	printf("-h:     Displays this help (Optional)\n");
	printf("-help:  Displays this help (Optional)\n");
	printf("--help: Displays this help (Optional)\n");
	printf("-?:     Displays this help (Optional)\n");
}

void ShowReCoDeServerUsage() {
	printf("Usage:\n");
	printf("ReCoDe -d ARG -p ARG -o ARG [-v ARG] [-?] [--help]\n");
	printf("\n");
	printf("-i:     (Required) Output image filename\n");
	printf("-o:     (Required) Output directory\n");
	printf("-d:     (Required) Dark file\n");
	printf("-p:     (Required) Params file\n");
	printf("-v:     Verbosity level (0 or 1, Optional)\n");
	printf("-l:     Log file name (Optional)\n");
	printf("-n:     Run name (Optional). Used while logging.\n");
	printf("-vf:    Gap between validation frames. (Optional). If not specified no validation frames are saved.\n");
	printf("-h:     Displays this help (Optional)\n");
	printf("-help:  Displays this help (Optional)\n");
	printf("--help: Displays this help (Optional)\n");
	printf("-?:     Displays this help (Optional)\n");
}

static const TCHAR *GetLastErrorText(int a_nError) {
	switch (a_nError) {
		case SO_SUCCESS:          return _T("Success");
		case SO_OPT_INVALID:      return _T("Unrecognized option");
		case SO_OPT_MULTIPLE:     return _T("Option matched multiple strings");
		case SO_ARG_INVALID:      return _T("Option does not accept argument");
		case SO_ARG_INVALID_TYPE: return _T("Invalid argument format");
		case SO_ARG_MISSING:      return _T("Required argument is missing");
		case SO_ARG_INVALID_DATA: return _T("Invalid argument data");
		default:                  return _T("Unknown error");
	}
}

int parse_recode_server_init_params (int argc, TCHAR * argv[], InitParams **init_params) {
	
	// declare our options parser, pass in the arguments from main
	// as well as our array of valid options.
	CSimpleOpt args(argc, argv, g_rgOptions);

	(*init_params)->verbosity = 0;
	(*init_params)->log_filename = "recode.log";
	(*init_params)->run_name = "Run";
	(*init_params)->validation_frame_gap = -1;

	// while there are arguments left to process
	while (args.Next()) {

		if (args.LastError() != SO_SUCCESS) {
			_tprintf(_T("%s: '%s' (use --help to get command line help)\n"), GetLastErrorText(args.LastError()), args.OptionText());
			continue;
		}

		switch (args.OptionId()) {
			case OPT_HELP:
				ShowReCoDeServerUsage();
				return 0;
			case ARG_VERB:
				(*init_params)->verbosity = wchar_to_int(args.OptionArg());
				break;
			case ARG_IMG:
				(*init_params)->image_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_DRK:
				(*init_params)->dark_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_PRM:
				(*init_params)->params_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_DIR:
				(*init_params)->output_directory = wchar_to_char(args.OptionArg());
				break;
			case ARG_LOG:
				(*init_params)->log_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_NAM:
				(*init_params)->run_name = wchar_to_char(args.OptionArg());
				break;
			case ARG_VAF:
				(*init_params)->validation_frame_gap = wchar_to_int(args.OptionArg());
				break;
			default:
				return 0;
		}
	}

	int params_check_cleared = 1;
	if ((*init_params)->output_directory == NULL) {
		printf("Missing Parameter Output Directory: Use -o option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->image_filename == NULL) {
		printf("Missing Parameter Image File: Use -i option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->mode == 0 && (*init_params)->dark_filename == NULL) {
		printf("Missing Parameter Dark File: Use -o option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->mode == 0 && (*init_params)->params_filename == NULL) {
		printf("Missing Parameter Params File: Use -p option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->verbosity > 1) {
		(*init_params)->verbosity = 1;
	}

	if (params_check_cleared == 0) {
		ShowReCoDeServerUsage();
		exit(0);
	}
	
	recode_print("ReCoDe Server Input Parameters:\n");
	recode_print("image_filename = %s\n", (*init_params)->image_filename);
	recode_print("verbosity = %d\n", (*init_params)->verbosity);
	recode_print("output_directory = %s\n", (*init_params)->output_directory);
	recode_print("dark_filename = %s\n", (*init_params)->dark_filename);
	recode_print("params_filename = %s\n", (*init_params)->params_filename);
	recode_print("log_filename = %s\n", (*init_params)->log_filename);
	recode_print("run_name = %s\n", (*init_params)->run_name);
	recode_print("validation_frame_gap = %d\n", (*init_params)->validation_frame_gap);
}

int parse_init_params(int argc, TCHAR * argv[], InitParams **init_params) {

	// declare our options parser, pass in the arguments from main
	// as well as our array of valid options.
	CSimpleOpt args(argc, argv, g_rgOptions);

	(*init_params)->verbosity = 0;
	(*init_params)->log_filename = "recode.log";
	(*init_params)->run_name = "Run";
	(*init_params)->validation_frame_gap = -1;
	(*init_params)->mode = -1;

	// while there are arguments left to process
	while (args.Next()) {

		if (args.LastError() != SO_SUCCESS) {
			_tprintf(_T("%s: '%s' (use --help to get command line help)\n"), GetLastErrorText(args.LastError()), args.OptionText());
			continue;
		}

		switch (args.OptionId()) {
			printf("%s", args.OptionId());
			case OPT_HELP:
				ShowUsage();
				return 0;
			case OPT_RC:
				if ((*init_params)->mode == -1) {
					(*init_params)->mode = 0;
				}
				else {
					printf("-rc and -de options cannot be used simultaneously.\n");
					return 0;
				}
				break;
			case OPT_DE:
				if ((*init_params)->mode == -1) {
					(*init_params)->mode = 1;
				}
				else {
					printf("-rc and -de options cannot be used simultaneously.\n");
					return 0;
				}
				break;
			case ARG_VERB:
				(*init_params)->verbosity = wchar_to_int(args.OptionArg());
				break;
			case ARG_IMG:
				(*init_params)->image_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_DRK:
				(*init_params)->dark_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_PRM:
				(*init_params)->params_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_DIR:
				(*init_params)->output_directory = wchar_to_char(args.OptionArg());
				break;
			case ARG_LOG:
				(*init_params)->log_filename = wchar_to_char(args.OptionArg());
				break;
			case ARG_NAM:
				(*init_params)->run_name = wchar_to_char(args.OptionArg());
				break;
			case ARG_VAF:
				(*init_params)->validation_frame_gap = wchar_to_int(args.OptionArg());
				break;
			default:
				return 0;
		}
	}

	int params_check_cleared = 1;
	if ((*init_params)->mode == -1) {
		printf("Missing Parameter: Specify -rc or -de.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->output_directory == NULL) {
		printf("Missing Parameter Output Directory: Use -o option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->image_filename == NULL) {
		printf("Missing Parameter Image File: Use -i option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->mode == 0 && (*init_params)->dark_filename == NULL) {
		printf("Missing Parameter Dark File: Use -o option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->mode == 0 && (*init_params)->params_filename == NULL) {
		printf("Missing Parameter Params File: Use -p option to specify.\n");
		params_check_cleared = 0;
	}
	if ((*init_params)->verbosity > 1) {
		(*init_params)->verbosity = 1;
	}

	if (params_check_cleared == 0) {
		ShowReCoDeServerUsage();
		exit(0);
	}

	/*
	recode_print("ReCoDe Server Input Parameters:\n");
	recode_print("image_filename = %s\n", (*init_params)->image_filename);
	recode_print("verbosity = %d\n", (*init_params)->verbosity);
	recode_print("output_directory = %s\n", (*init_params)->output_directory);
	recode_print("dark_filename = %s\n", (*init_params)->dark_filename);
	recode_print("params_filename = %s\n", (*init_params)->params_filename);
	recode_print("log_filename = %s\n", (*init_params)->log_filename);
	recode_print("run_name = %s\n", (*init_params)->run_name);
	recode_print("validation_frame_gap = %d\n", (*init_params)->validation_frame_gap);
	*/
}



