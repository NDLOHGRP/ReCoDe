void logResult (const char* testname, float* compress_times_sizes) {
	
	FILE *logfile = fopen("results.txt", "a");
	float compression_ratio = compress_times_sizes[3]/compress_times_sizes[2];
	fprintf(logfile, "%s %2f %2f %2f %2f %2f\n", testname, compress_times_sizes[0], compress_times_sizes[1], compress_times_sizes[2], compress_times_sizes[3], compression_ratio);
	fclose(logfile);
}

void logResultHeader (const char* expname, const char *fieldnames[], int nFields) {

	FILE *logfile1 = fopen("results.txt", "a");
	if (logfile1 == NULL)
	{
		printf("Error opening log file!\n");
	}
	fprintf(logfile1, "\n");
	
	fprintf(logfile1, "Experiment Name : %s\n", expname);
	int i;
	for (i=0; i<nFields; i++) {
		fprintf(logfile1, "%s ", fieldnames[i]);
	}
	
	fprintf(logfile1, "\n");
	fclose(logfile1);
	
}