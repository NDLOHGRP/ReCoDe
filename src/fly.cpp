/**************************************************************************************************************
to compile:    gcc fly.cpp -L../zlib-1.2.11 -lz -L../bzip2-1.0.6 -lbz2 -lsnappy -L../lz4-1.8.2/lib -llz4 -O3 -o fly -lm -fopenmp
to compile with profiling: kinst-ompp-papi gcc multifly.c -lz -O3 -o fly -lm -fopenmp

to run:        ./fly

-lz:     links zlib
-lbz2:     links bzip2
-lm:     links math.h NOTE: it is important to keep -lm as the last option
-o:         write build output to file
-O3:     optimizer flag
**************************************************************************************************************/

#ifndef ENABLE_MULTIPLE_COMPRESSIONS
#define ENABLE_MULTIPLE_COMPRESSIONS 0
#endif

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <zmq.h>
#include <assert.h>

#include "zlib.h"

#if ENABLE_MULTIPLE_COMPRESSIONS
#include "bzlib.h"
#include "lzma.h"
#include "snappy-c.h"
#include "lz4.h"
#endif

#include <tchar.h>
#include <wchar.h>
#include <SimpleOpt.h>

#ifdef _WIN32
#include <win32\dirent.h>
#include <windows.h>
#else
#include <dirent.h>
#endif

#include "recodefs.h"
#include "initparser.h"
#include "recode_utils.h"
#include "mrchandler.h"
#include "sequence_reader.h"
#include "commons.h"
#include "compressor.h"
#include "logger.h"
#include "fileutils.h"
/*
#include "mrc_header.h"
*/
#include "argparser.h"
#include "recode_header.h"
/*
#include "ReCoDeConfig.h"
*/
#include "L1.h"
#include "L2.h"
#include "L3.h"
#include "L4.h"

#define    RC_WRITE_MODE_NO_RETURN      0        // reduceCompressFrame_L1 but do not write to return buffer
#define    RC_WRITE_MODE_RETURN_ONLY	1        // reduceCompressFrame_L1 and writes to return buffer but the data is not copied to OLT buffer
#define    RC_WRITE_MODE_RETURN_COPY	2        // reduceCompressFrame_L1 and writes to return buffer and the data is copied to OLT buffer

// Global Variables
unsigned int	WAIT_MILLISECONDS	= 1;
char*			IMAGE_FILENAME		= "";
RunStats*		RUN_STATS;
uint8_t*		RCT_STATE_INDICATORS;
uint8_t			RUN_STATE_INDICATOR = 1;		// 1 = RUNNING; 2 = SHUTTING DOWN

// reduce-compress thread
void runRCT_L1(
	InitParams *init_params,
	InputParams *params,
	RCHeader *header,
	uint16_t* darkFrame,
	uint8_t process_id,
	uint8_t MPI_PID_Offset,
	uint8_t *pBuffer,
	uint64_t BUFFER_LENGTH
) {

	printf("RCT: %d: Starting... \n", process_id);

	uint8_t rc_operation_mode = params->reduction_level;
	uint8_t rc_write_mode = 2;
	char* destination_directory = init_params->output_directory;
	char* name = init_params->run_name;
	int validation_frame_gap = init_params->validation_frame_gap;
	uint8_t num_threads_in_pool = params->num_threads;

	uint64_t MAX_FRAMES_PER_CHUNK = params->num_frames;
	uint64_t MAX_FRAMES_PER_THREAD = ceil((MAX_FRAMES_PER_CHUNK*1.0) / (params->num_threads*1.0));

	uint64_t i;
	uint64_t run_count = 0;
	uint64_t count = 0;
	uint64_t wait_time = 0;
	long available_buffer_space = BUFFER_LENGTH;
	uint64_t RCT_Buffer_Fill_Position = 0;
	uint8_t RCT_FULL_IND = 0;
	uint64_t nWrites = 0;
	uint64_t run_offset = 0;
	uint64_t absolute_frame_index;
	uint64_t n_frames_in_chunk;

	uint8_t copy_compressed_frame_to_return_buffer;
	if (rc_write_mode > RC_WRITE_MODE_NO_RETURN) {
		copy_compressed_frame_to_return_buffer = 1;
	} else {
		copy_compressed_frame_to_return_buffer = 0;
	}

	uint32_t n_bytes_in_packed_pixvals, n_compressed_bytes_1, n_compressed_bytes_2, compressed_frame_length;

	// create data buffers
	uint32_t n_pixels_in_frame = header->nx * header->ny;
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);
	uint16_t *pixvals = (uint16_t*)malloc(n_pixels_in_frame * sizeof(uint16_t));        // allocated max space - where all pixels are fg pixels, must calloc to init to 0
	uint8_t  *binaryImage = (uint8_t*)malloc(n_bytes_in_binary_image * sizeof(uint8_t));    // every bit in binaryImage has to be reset for every frame in the loop below
	uint8_t  *packedPixvals = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressedPixvals = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressedBinaryImage = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressedPackedPixvals = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressed_frame = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));

	const char* filename = makePartFilename(process_id, name, 1);
	char* bfilename = concat(destination_directory, filename);
	printf("RCT: %d: Will be using part file '%s'\n", process_id, bfilename);

	FILE *partFile = fopen(bfilename, "wb");
	if (!partFile) {
		printf("Failed to create file: %s.", bfilename);
		exit(1);
	}
	int FD = fileno(partFile);
	fwrite(&process_id, sizeof(char), 1, partFile);        // write process_id to partfile

	FILE *validationFile;
	if (validation_frame_gap > 0) {
		char* vfilename = concat(concat(destination_directory, makeValidationFilename(process_id, name)), ".bin");
		validationFile = fopen(vfilename, "wb");
		if (!validationFile) {
			printf("Failed to create file: %s.", vfilename);
			exit(1);
		}
	}

	float *run_metrics = (float*)malloc(5 * sizeof(float));
	uint64_t sz_frameBuffer = header->ny * header->nx * MAX_FRAMES_PER_THREAD * sizeof(uint16_t);
	uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);

	while (1) {

		//printf("RCT %d: is alive\n", process_id);

		if (RCT_STATE_INDICATORS[process_id] == 0) {
			continue;
		}

		/*
		if (RCT_STATE_INDICATORS[process_id] == 1) {
			Sleep(1000);
			RCT_STATE_INDICATORS[process_id] = 0;
			run_count++;
			continue;
		}

		if (RUN_STATE_INDICATOR == 0) {
			return;
		}
		*/

		if (RCT_STATE_INDICATORS[process_id] == 1) {

			// Load header for this thread
			SEQHeader *seq_header = (SEQHeader *)malloc(sizeof(SEQHeader));
			FILE* fp = fopen(IMAGE_FILENAME, "rb");
			if (!fp) {
				printf("Failed to open file: '%s' for reading.", IMAGE_FILENAME);
				exit(1);
			}
			getSEQHeader(fp, &seq_header);

			// Determine the number of available frames for this thread
			n_frames_in_chunk = seq_header->AllocatedFrames;
			uint64_t n_frames_per_thread = ceil((n_frames_in_chunk*1.0) / (params->num_threads*1.0));
			uint64_t frame_offset = process_id * n_frames_per_thread;
			uint64_t num_frames_to_process = min(n_frames_per_thread, n_frames_in_chunk - frame_offset);

			// Load data for this thread
			uint64_t available_frames = getSEQFrames(process_id, fp, frameBuffer, frame_offset, num_frames_to_process, seq_header);
			fclose(fp);

			recode_print("RCT %d: Available Frames: %d; Num. Frames to Process: %d\n", process_id, available_frames, num_frames_to_process);

			// process frames
			while (count < available_frames) {

				absolute_frame_index = run_offset + frame_offset + count;
				if (validation_frame_gap > 0) {
					if (absolute_frame_index%validation_frame_gap == 0) {
						fwrite(&absolute_frame_index, sizeof(uint64_t), 1, validationFile);
						fwrite(&frameBuffer[count*n_pixels_in_frame], sizeof(uint16_t), n_pixels_in_frame, validationFile);
						recode_print("RCT %d: Saving frame: %d\n", process_id, absolute_frame_index);
					}
				}
				//recode_print("RCT %d: Processing frame: %d of %d [frame index = %d]\n", process_id, count, num_frames_to_process, absolute_frame_index);
				count++;

				/*
				reduceCompressFrame_L1(
					process_id,
					frameBuffer,
					darkFrame,
					absolute_frame_index,                // absolute frame index
					count,								// frame index relative to start of this thread
					h,
					rc_operation_mode,                    // rc_operation_mode: 1 = reduce and compress
					params->dark_threshold_epsilon,                                    // epsilon_s
					params->source_bit_depth,                                    // bit_depth
					0,                                    // compression_scheme: 0 = gzip
					1,                                    // compression_level: 1 = fastest, 9 = best compression
					pixvals,
					binaryImage,
					packedPixvals,
					compressedBinaryImage,
					compressedPackedPixvals,
					&n_bytes_in_packed_pixvals,
					&n_compressed_bytes_1,
					&n_compressed_bytes_2,
					run_metrics,
					copy_compressed_frame_to_return_buffer,    // copy_compressed_frame_to_return_buffer
					compressed_frame,
					&compressed_frame_length
				);

				if (available_buffer_space >= compressed_frame_length) {
					if (rc_write_mode > RC_WRITE_MODE_RETURN_ONLY) {
						for (i = 0; i < compressed_frame_length; i++) {
							pBuffer[RCT_Buffer_Fill_Position + i] = compressed_frame[i];
						}
						//recode_print("Copying to OLT Buffer.\n");
					}

					RCT_Buffer_Fill_Position += compressed_frame_length;
					if (RCT_Buffer_Fill_Position > BUFFER_LENGTH) {
						printf("RCT %d: Buffer Fill Position = %llu", process_id, RCT_Buffer_Fill_Position);
					}
					available_buffer_space -= compressed_frame_length;
					recode_print("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
					count++;

				}
				else {

					recode_print("RCT %d: buffer is full, waiting for OLT to clear.\n", process_id);
					fwrite(pBuffer, sizeof(char), RCT_Buffer_Fill_Position, partFile);
					nWrites++;

					available_buffer_space = BUFFER_LENGTH;
					RCT_Buffer_Fill_Position = 0;
					RCT_FULL_IND = 0;
					recode_print("RCT %d: OLT cleared buffer.\n", process_id);

					// copy the last received frame that was never written
					if (rc_write_mode > RC_WRITE_MODE_RETURN_ONLY) {
						for (i = 0; i < compressed_frame_length; i++) {
							pBuffer[RCT_Buffer_Fill_Position + i] = compressed_frame[i];
						}
						//recode_print("Copying to OLT Buffer.\n");
						recode_print("RCT %d \t Waited Milliseconds: \t%d.\n", process_id, WAIT_MILLISECONDS*wait_time);
						wait_time = 0;
					}

					RCT_Buffer_Fill_Position += compressed_frame_length;
					if (RCT_Buffer_Fill_Position > BUFFER_LENGTH) {
						printf("RCT %d: Buffer Fill Position = %llu", process_id, RCT_Buffer_Fill_Position);
					}
					available_buffer_space -= compressed_frame_length;
					count++;
					recode_print("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
				}
				*/
			}

			count = 0;
			run_offset += n_frames_in_chunk;
			RCT_STATE_INDICATORS[process_id] = 0;

		}

		if (RUN_STATE_INDICATOR == 0) {

			RCT_STATE_INDICATORS[process_id] = 1;

			RCT_FULL_IND = 3;                            // this RCT is done, indicate OLT to close RCT's file

			fclose(partFile);
			if (validation_frame_gap > 0) {
				fclose(validationFile);
			}

			free(pixvals);
			free(binaryImage);
			free(packedPixvals);
			free(compressedPixvals);
			free(compressedBinaryImage);
			free(compressedPackedPixvals);
			free(compressed_frame);
			recode_print("RCT %d: Done (total wait time = %u milliseconds).\n", process_id, wait_time * WAIT_MILLISECONDS);

			RUN_STATS[process_id].thread_id = process_id;
			RUN_STATS[process_id].reduction_time = run_metrics[0];
			RUN_STATS[process_id].compression_time = run_metrics[1];
			RUN_STATS[process_id].reduced_bytes = run_metrics[2];
			RUN_STATS[process_id].compressed_bytes = run_metrics[3];
			RUN_STATS[process_id].decompression_time = run_metrics[4];
			RUN_STATS[process_id].num_writes = nWrites;
			RUN_STATS[process_id].wait_time = wait_time * WAIT_MILLISECONDS;

			RCT_STATE_INDICATORS[process_id] = 0;

			return;
		}

	}
	
}

void start_recode_server(InitParams *init_params, InputParams *params, RCHeader *header, uint16_t* darkFrame) {

	uint8_t MPI_PID = 0;

	int i = 0;

	/*
	==========================================================================================
	Create RCT buffers
	==========================================================================================
	*/
	uint64_t BUFFER_LENGTH = params->num_frames*params->num_rows*params->num_cols;
	uint8_t **pBuffers;
	pBuffers = (uint8_t**)malloc(params->num_threads * sizeof(uint8_t*));
	for (i = 0; i < params->num_threads; i++) {
		pBuffers[i] = (uint8_t *)malloc(BUFFER_LENGTH * sizeof(uint8_t)); // make actual arrays
	}

	/*
	==========================================================================================
	Do multithreaded reduction-compression
	==========================================================================================
	*/

	/*
	time_t timer;
	char buffer[27];
	struct tm* tm_info;
	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer, 27, "%d-%b-%Y %H:%M:%S", tm_info);
	*/

	RUN_STATS = (RunStats*)calloc(params->num_threads, sizeof(RunStats));
	for (i = 0; i < params->num_threads; i++) {
		RUN_STATS[i].run_name = init_params->run_name;
		RUN_STATS[i].start_time = "";
	}

	#pragma omp parallel for num_threads(params->num_threads)
	for (i = 0; i < params->num_threads; i++) {
		if (params->reduction_level == 1) {
			runRCT_L1(
				init_params,
				params,
				header,
				darkFrame,
				i,
				MPI_PID*params->num_threads,
				pBuffers[i],
				BUFFER_LENGTH
			);
		}
	}

}

int get_RCT_States(uint8_t *rct_state_indicators, uint8_t nThreads, uint8_t val) {
	int i;
	int count = 0;
	for (i = 0; i < nThreads; i++) {
		if (rct_state_indicators[i] == val) {
			count++;
		}
	}
	return count;
}

void set_RCT_States(uint8_t *rct_state_indicators, uint8_t nThreads, uint8_t val) {
	int i;
	for (i = 0; i < nThreads; i++) {
		rct_state_indicators[i] = val;
	}
}

void log_last_run_stats(char* log_filename, InputParams *params) {

}

void reset_run_stats() {

}

int start_zmq_server(InitParams *init_params, InputParams *params) {

	//  Socket to talk to clients
	void *context = zmq_ctx_new();
	void *responder = zmq_socket(context, ZMQ_REP);
	int rc = zmq_bind(responder, "tcp://127.0.0.1:18534");
	assert(rc == 0);

	int request_count = 0;

	RUN_STATE_INDICATOR = 1;

	while (1) {
		char buffer[500];
		zmq_recv(responder, buffer, 500, 0);
		printf("ReCoDe Server: Received Request %d, processing: %s\n", request_count, buffer);
		log_last_run_stats(init_params->log_filename, params);
		if (strncmp(buffer, "shutdown", 8) == 0) {
			RUN_STATE_INDICATOR = 0;
			Sleep(100);          //  Wait for RCT Threads to finish
			while (get_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1) != 0) {
				Sleep(100);          //  Wait for RCT Threads to finish
			}
			zmq_send(responder, "Shutdown Complete", 17, 0);
			exit(0);
		}
		reset_run_stats();
		IMAGE_FILENAME = buffer;
		set_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1);
		Sleep(100);
		while (get_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1) != 0) {
			Sleep(100);          //  Wait for RCT Threads to finish
		}
		zmq_send(responder, "Request Completed", 17, 0);
		request_count++;
	}
	return 0;
}

void _simulate_zmq_run(InitParams *init_params, InputParams *params) {

	int iters = 3;
	for (int i = 0; i < iters; i++) {
		printf("ReCoDe Server: Received Request %d, processing: %d\n", i, i);
		reset_run_stats();
		IMAGE_FILENAME = "R:\\003.seq";
		set_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1);
		Sleep(3000);
		while (get_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1) != 0) {
			Sleep(100);          //  Wait for RCT Threads to finish
		}
		printf("Request Completed\n");
	}
	
	RUN_STATE_INDICATOR = 0;
	while (get_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1) != 0) {
		Sleep(100);          //  Wait for RCT Threads to finish
	}
	printf("Shutdown Complete\n");
	exit(0);

}

#ifdef _WIN32
int _tmain(int argc, TCHAR * argv[]) {

	/*
	const char* filename = "D:/cbis/images/20161207/14-37-50.811.seq";
	DataSize d = { 4096, 512, 50, 12 };
	uint64_t sz_darkBuffer = d.nx * d.ny * d.nz * sizeof(uint16_t);
	uint16_t *darkBuffer = (uint16_t *)malloc(sz_darkBuffer);
	
	SEQHeader *seq_header = (SEQHeader *)malloc(sizeof(SEQHeader));
	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		printf("Failed to open file: '%s' for reading.", filename);
		exit(1);
	}
	getSEQHeader(fp, &seq_header);
	uint64_t available_frames = getSEQFrames(fp, darkBuffer, 300, 50, seq_header);
	fclose(fp);
	
	createSEQFiles(darkBuffer, seq_header, 50, "D:/cbis/GitHub/ReCoDe/scratch/temp/", 10);
	return 0;
	*/

	omp_set_nested(1);

	InitParams *init_params = (InitParams *)malloc(sizeof(InitParams));
	parse_recode_server_init_params(argc, argv, &init_params);

	// read and validate input parameters
	InputParams *params = (InputParams *)malloc(sizeof(InputParams));
	get_input_params(init_params->params_filename, &params);

	// if MRC file, replace params from MRC header
	if (params->source_file_type == 1) {
		MRCHeader *mrc_header = (MRCHeader *)malloc(sizeof(MRCHeader));
		parseMRCHeader(init_params->image_filename, &mrc_header);
		compile_missing_params(mrc_header, &params);
		if (params->source_bit_depth != 16) {
			printf("ReCoDe currently only supports 16-bit unsigned data. The provided data has bit-depth %d.", params->source_bit_depth);
			exit(1);
		}
		free(mrc_header);
	}

	/*
	==========================================================================================
	Create ReCoDE header
	==========================================================================================
	*/
	RCHeader *header = (RCHeader *)malloc(sizeof(RCHeader));
	uint8_t *image_name = (uint8_t*)getFilenameFromPath(init_params->image_filename);
	uint8_t *dark_name = (uint8_t*)getFilenameFromPath(init_params->dark_filename);
	create_recode_header(params, -1, image_name, dark_name, &header);
	print_recode_header(header);

	/*
	==========================================================================================
	Dark Noise Estimation
	==========================================================================================
	*/
	uint32_t n_pixels_in_frame = header->nx * header->ny;
	uint16_t *darkFrame = (uint16_t *)calloc(n_pixels_in_frame, sizeof(uint16_t));

	recode_print("header->num_dark_frames = %lu\n", header->num_dark_frames);

	if (header->num_dark_frames > 1) {

		DataSize d = { header->nx, header->ny, params->num_dark_frames, params->source_bit_depth };
		uint64_t sz_darkBuffer = d.nx * d.ny * d.nz * sizeof(uint16_t);
		uint16_t *darkBuffer = (uint16_t *)malloc(sz_darkBuffer);

		loadData(params->dark_file_type, init_params->dark_filename, darkBuffer, 0, header->num_dark_frames, d);
		//serializeFrames(darkBuffer, header->nx, header->ny, "dark_frame.txt", 1);
		printf("Dark data loaded.\n");

		getDarkMax(darkBuffer, d, darkFrame, params->num_threads);
		free(darkBuffer);

		printf("Dark noise estimation complete.\n");

	}
	else {
		
		DataSize d1 = { header->nx, header->ny, 1, params->source_bit_depth };
		loadData(params->dark_file_type, init_params->dark_filename, darkFrame, 0, 1, d1);
		//serializeFrames(darkFrame, header->nx, header->ny, "dark_frame.txt", 1);
		printf("Pre-computed dark noise estimates loaded.\n");

	}

	/* RCT State Indicators: shared between RCTs and ZMQ */
	RCT_STATE_INDICATORS = (uint8_t*)calloc(params->num_threads, sizeof(uint8_t));
	for (int i = 0; i<params->num_threads; i++) {
		RCT_STATE_INDICATORS[i] = 0;
	}

	#pragma omp parallel for num_threads(2)
	for (int i = 0; i < 2; i++) {
		if (i == 0) {
			start_recode_server (init_params, params, header, darkFrame);
		}
		else {
			start_zmq_server (init_params, params);
			//_simulate_zmq_run(init_params, params);
		}
	}

	/*
	==========================================================================================
	Cleanup
	==========================================================================================
	*/
	recode_print("0\n");
	free(darkFrame);
	recode_print("1\n");
}
#else
int main(int argc, char *argv[]) {

	CSimpleOpt args(argc, argv, g_rgOptions);

	char* name = argv[1];
	uint64_t buffer_size = atoi(argv[2]);
	uint64_t nImageFrames = atoi(argv[3]);
	uint16_t min_threads = atoi(argv[4]);
	uint16_t max_threads = atoi(argv[5]);
	uint16_t threads_step = atoi(argv[6]);
	uint8_t dose_rate = atof(argv[7]);
	uint8_t rc_operation_mode = atoi(argv[8]);	// rc_operation_mode: 1 = reduce and compress; 0 = reduce only
	uint8_t replicates = atoi(argv[9]);
	char* imageFile = argv[10];
	char* darkNoiseFile = argv[11];
	char* destination_directory = argv[12];

	DataSize b = { 4096, 512, nImageFrames, 16 };
	uint16_t *frameBuffer = (uint16_t *)malloc(b.nx * b.ny * b.nz * sizeof(uint16_t));
	loadData(imageFile, frameBuffer, 0, nImageFrames, b);
	recode_print("Image data loaded.\n");

	/*
	==========================================================================================
	Dark Noise Estimation
	==========================================================================================
	*/
	uint16_t *darkFrame = (uint16_t *)calloc(4096 * 512, sizeof(uint16_t));
	DataSize d1 = { 4096, 512, 1, 12 };
	loadData(darkNoiseFile, darkFrame, 0, 1, d1);
	recode_print("Pre-computed dark noise estimates loaded.\n");

	/*
	==========================================================================================
	Create Log
	==========================================================================================
	*/
	time_t timer;
	char buffer[27];
	struct tm* tm_info;

	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer, 27, "%d-%b-%Y %H:%M:%S", tm_info);

	FILE *logfile = fopen("recode_on_the_fly_results", "a");

	fprintf(logfile, "\nRun Start Time: %s\n", buffer);
	fprintf(logfile, "Data Saved to NAS by Each RC Thread independently (DE-16)\n");
	fprintf(logfile, "=======================================\n");
	fprintf(logfile, "L \t DR \t NT \t n_Fr \t B_Sz \t T \t R/C \t WMODE \t REP\n");
	fclose(logfile);

	uint8_t reps;
	uint8_t level = 1;
	uint8_t rc_write_mode = 3;
	uint8_t rct_threads;

	int taskid = -1;

	for (rct_threads = min_threads; rct_threads < max_threads + 1; rct_threads = rct_threads + threads_step) {
		for (reps = 0; reps < replicates; reps++) {

			recode_on_the_fly(
				level,
				rc_operation_mode,
				rc_write_mode,
				frameBuffer,
				darkFrame,
				rct_threads,
				buffer_size,
				nImageFrames,
				0,
				dose_rate,
				taskid,
				reps,
				destination_directory,
				name
			);

		}
	}


	/*
	==========================================================================================
	Cleanup
	==========================================================================================
	*/
	recode_print("0\n");
	free(darkFrame);
	free(frameBuffer);
	recode_print("1\n");
}
#endif

