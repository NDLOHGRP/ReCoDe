/**************************************************************************************************************
to compile:    gcc fly.cpp -L../zlib-1.2.11 -lz -L../bzip2-1.0.6 -lbz2 -lsnappy -L../lz4-1.8.2/lib -llz4 -O3 -o fly -lm -fopenmp
to compile with profiling: kinst-ompp-papi gcc multifly.c -lz -O3 -o fly -lm -fopenmp

to run:        ./fly

-lz:     links zlib
-lbz2:   links bzip2
-lm:     links math.h NOTE: it is important to keep -lm as the last option
-o:      write build output to file
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
// unsigned int	WAIT_MILLISECONDS	= 1;
char*			IMAGE_FILENAME		= "";
RunStats*		RUN_STATS;
uint8_t*		RCT_STATE_INDICATORS;
uint8_t			RUN_STATE_INDICATOR = 1;		// 1 = RUNNING; 2 = SHUTTING DOWN

// reduce-compress thread
void runRCT_L1(
	InitParams *init_params,
	InputParams *params,
	RCHeader *header,
	uint16_t *darkFrame,
	uint8_t process_id,
	uint8_t MPI_PID_Offset,
	uint8_t *pBuffer,
	uint64_t BUFFER_LENGTH
) {

	printf("RCT: %d: Starting... \n", process_id);

	uint8_t copy_compressed_frame_to_return_buffer = 1;
	uint8_t rc_operation_mode = params->reduction_level;
	char* destination_directory = init_params->output_directory;
	char* name = init_params->run_name;
	int validation_frame_gap = init_params->validation_frame_gap;
	uint8_t num_threads_in_pool = params->num_threads;

	uint64_t MAX_FRAMES_PER_CHUNK = params->num_frames;
	uint64_t MAX_FRAMES_PER_THREAD = ceil((MAX_FRAMES_PER_CHUNK*1.0) / (params->num_threads*1.0));

	uint64_t i;
	uint64_t run_count = 0;
	uint64_t count = 0;
	uint32_t num_frames_in_part = 0;
	uint64_t wait_time = 0;
	uint64_t available_buffer_space = BUFFER_LENGTH;
	uint64_t RCT_Buffer_Fill_Position = 0;
	uint8_t RCT_FULL_IND = 0;
	uint64_t nWrites = 0;
	uint64_t run_offset = 0;
	uint64_t absolute_frame_index, n_frames_in_chunk;
	uint32_t n_bytes_in_packed_pixvals, n_compressed_bytes_1, n_compressed_bytes_2, compressed_frame_length;

	// create data buffers
	uint32_t n_pixels_in_frame = header->nx * header->ny;
	uint32_t n_bytes_in_binary_image = ceil(n_pixels_in_frame / 8.0);
	uint16_t *pixvals = (uint16_t*)malloc(n_pixels_in_frame * sizeof(uint16_t));			// allocated max space - where all pixels are fg pixels, must calloc to init to 0
	uint8_t  *binaryImage = (uint8_t*)malloc(n_bytes_in_binary_image * sizeof(uint8_t));    // every bit in binaryImage has to be reset for every frame in the loop below
	uint8_t  *packedPixvals = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressedPixvals = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressedBinaryImage = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressedPackedPixvals = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));
	uint8_t  *compressed_frame = (uint8_t*)malloc(n_pixels_in_frame * sizeof(uint8_t));

	/*============== temporary variables for counting on validation frames ===========
	only used if params->reduction_level == 1 or params->reduction_level == 3
	all variables usd in validation counting use _vc prefix
	=================================================================================*/
	uint64_t _vc_frame_start_index;
	uint32_t _vc_n_fg_pixels, _vc_n_labels;
	DataSize _vc_h = { 128, 128, 1, params->source_bit_depth };
	if (header->nx < 128) {
		_vc_h.nx = header->nx;
	}
	if (header->ny < 128) {
		_vc_h.ny = header->ny;
	}
	uint32_t _vc_left = (uint32_t)floor((header->nx - _vc_h.nx) / 2.0);
	uint32_t _vc_top = (uint32_t)floor((header->ny - _vc_h.ny) / 2.0);
	uint32_t _vc_roi_top_left = header->nx * _vc_top + _vc_left;
	uint32_t _vc_n_pixels_in_frame = _vc_h.nx * _vc_h.ny;
	float *_vc_x_coor = (float *)malloc(_vc_n_pixels_in_frame * sizeof(float));
	float *_vc_y_coor = (float *)malloc(_vc_n_pixels_in_frame * sizeof(float));
	uint16_t *_vc_foregroundImage = (uint16_t *)calloc(_vc_n_pixels_in_frame, sizeof(uint16_t));
	int8_t   *_vc_foregroundTernaryMap = (int8_t *)calloc(_vc_n_pixels_in_frame, sizeof(int8_t));
	/*============== temporary variables for counting on validation frames ===========*/

	const char* filename = makePartFilename(process_id, name, params->reduction_level);
	char* bfilename = concat(destination_directory, filename);
	printf("RCT: %d: Will be using part file '%s'\n", process_id, bfilename);

	FILE *partFile = fopen(bfilename, "wb");
	if (!partFile) {
		printf("Failed to create file: %s.", bfilename);
		exit(1);
	}
	int FD = fileno(partFile);

	// serialize header
	serialize_recode_header(partFile, header);
	//fwrite(&process_id, sizeof(char), 1, partFile);        // write process_id to partfile

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

			if (RUN_STATE_INDICATOR == 0) {

				RCT_STATE_INDICATORS[process_id] = 1;

				recode_print("RCT %d: Clearing buffer and closing files...\n", process_id);
				fwrite(pBuffer, sizeof(char), RCT_Buffer_Fill_Position, partFile);
				//fwrite(&num_frames_in_part, sizeof(uint32_t), 1, partFile);		// serialize the number of frames at the end
				nWrites++;

				// serialize the number of frames in the header
				fseek(partFile, 17, SEEK_SET);
				fwrite(&num_frames_in_part, sizeof(uint32_t), 1, partFile);
				fseek(partFile, 277, SEEK_SET);
				fwrite(&process_id, sizeof(uint8_t), 1, partFile);

				fclose(partFile);
				if (validation_frame_gap > 0) {
					fclose(validationFile);
				}
				recode_print("RCT %d: Done.\n", process_id);

				recode_print("RCT %d: Freeing memory...\n", process_id);
				free(pixvals);
				free(binaryImage);
				free(packedPixvals);
				free(compressedPixvals);
				free(compressedBinaryImage);
				free(compressedPackedPixvals);
				free(compressed_frame);

				RUN_STATS[process_id].thread_id = process_id;
				RUN_STATS[process_id].reduction_time = run_metrics[0];
				RUN_STATS[process_id].compression_time = run_metrics[1];
				RUN_STATS[process_id].reduced_bytes = run_metrics[2];
				RUN_STATS[process_id].compressed_bytes = run_metrics[3];
				RUN_STATS[process_id].decompression_time = run_metrics[4];
				RUN_STATS[process_id].num_writes = nWrites;
				RUN_STATS[process_id].wait_time = 0;

				RCT_STATE_INDICATORS[process_id] = 0;
				recode_print("RCT %d: Done.\n", process_id);

				return;
			}

			continue;
		}

		if (RCT_STATE_INDICATORS[process_id] == 1) {

			// Load header for this thread
			SEQHeader *seq_header = (SEQHeader *)malloc(sizeof(SEQHeader));
			FILE* fp = fopen(IMAGE_FILENAME, "rb");
			if (!fp) {
				printf("Failed to open file: '%s' for reading.", IMAGE_FILENAME);
				exit(1);
			}
			getSEQHeader(fp, &seq_header);

			// Do basic sanity checks
			uint64_t frame_sz_in_bytes = seq_header->ImageInfo.ImageWidth * seq_header->ImageInfo.ImageHeight * (seq_header->ImageInfo.ImageBitDepth / 8);
			if (header->nx != seq_header->ImageInfo.ImageWidth) {
				printf("Input parameter 'nx' (%d) does not match actual image width (%d). Quitting.\n", header->nx, seq_header->ImageInfo.ImageWidth);
				exit(1);
			}
			if (header->ny != seq_header->ImageInfo.ImageHeight) {
				printf("Input parameter 'ny' (%d) does not match actual image height (%d). Quitting.\n", header->nx, seq_header->ImageInfo.ImageHeight);
				exit(1);
			}
			if (params->source_bit_depth != seq_header->ImageInfo.ImageBitDepthReal) {
				printf("Input parameter 'bit_depth' (%d) does not match actual image bit depth (%d). Quitting.\n", header->bit_depth, seq_header->ImageInfo.ImageBitDepthReal);
				exit(1);
			}

			// Determine the number of available frames for this thread
			n_frames_in_chunk = seq_header->AllocatedFrames;
			uint64_t n_frames_per_thread = ceil((n_frames_in_chunk*1.0) / (params->num_threads*1.0));
			uint64_t frame_offset = process_id * n_frames_per_thread;
			uint64_t num_frames_to_process = min(n_frames_per_thread, n_frames_in_chunk - frame_offset);

			// Load data for this thread
			uint64_t available_frames = getSEQFrames(process_id, fp, frameBuffer, frame_offset, num_frames_to_process, seq_header);
			/*===============Testing===============
			uint64_t i, j, frame_start_index, linear_index;
			uint16_t max_val = 0;
			int x[] = { 99, 1111, 2222, 3333, 4000 };
			int y[] = { 110, 220, 330, 440, 500 };
			for (i = 0; i < available_frames; i++) {
				frame_start_index = i* header->nx * header->ny;
				printf("Frame %" PRIu64 ": ", i + frame_offset);
				for (j = 0; j < 5; j++) {
					linear_index = y[j] * header->nx + x[j];
					printf("%d,%d=%" PRIu16 ", ", y[j], x[j], frameBuffer[frame_start_index + linear_index]);
				}
				printf("\n");
			}
			/*===============Testing===============*/
			fclose(fp);

			recode_print("RCT %d: Available Frames: %d; Num. Frames to Process: %d\n", process_id, available_frames, num_frames_to_process);

			// process frames
			while (count < available_frames) {

				absolute_frame_index = run_offset + frame_offset + count;

				// serialize validation data
				if (validation_frame_gap > 0) {
					if (absolute_frame_index%validation_frame_gap == 0) {
						fwrite(&absolute_frame_index, sizeof(uint64_t), 1, validationFile);
						fwrite(&frameBuffer[count*n_pixels_in_frame], sizeof(uint16_t), n_pixels_in_frame, validationFile);
						recode_print("RCT %d: Saving frame: %d\n", process_id, absolute_frame_index);
						if (params->reduction_level == 1 || params->reduction_level == 3) {
							// Dark Subtraction
							_vc_n_fg_pixels = 0;
							_vc_frame_start_index = count*n_pixels_in_frame + _vc_roi_top_left;
							get_foreground_image(frameBuffer + _vc_frame_start_index, darkFrame + _vc_roi_top_left, params->dark_threshold_epsilon, _vc_h, &_vc_n_fg_pixels, _vc_foregroundImage, _vc_foregroundTernaryMap);
							// Connected Components Analysis to get (binary) Centroid Image
							uint32_t _temp_n_labels;
							cca(_vc_foregroundImage, _vc_foregroundTernaryMap, frameBuffer + _vc_frame_start_index, darkFrame + _vc_roi_top_left, params->dark_threshold_epsilon, _vc_h, _vc_n_fg_pixels, _vc_x_coor, _vc_y_coor, &_vc_n_labels);
							float _vc_dose_rate = (float)_vc_n_labels / (_vc_h.nx*_vc_h.ny);
							recode_print("RCT %d: Num. Electron Events in Selected RoI of Frame %d: %" PRIu32 ", Est. Dose Rate: %f e/pixel/frame\n", process_id, absolute_frame_index, _vc_n_labels, _vc_dose_rate);
						}
					}
				}
				
				// reduce compress
				DataSize h = { header->nx, header->ny, num_frames_to_process, header->bit_depth };
				reduceCompressFrame_L1(
					process_id,
					frameBuffer,
					darkFrame,
					absolute_frame_index,				// absolute frame index
					count,								// frame index relative to start of this thread
					h,
					rc_operation_mode,                  // rc_operation_mode: 1 = reduce and compress
					params->dark_threshold_epsilon,     // epsilon_s
					params->source_bit_depth,           // bit_depth
					0,                                  // compression_scheme: 0 = gzip
					1,                                  // compression_level: 1 = fastest, 9 = best compression
					pixvals,
					binaryImage,
					packedPixvals,
					compressedBinaryImage,
					compressedPackedPixvals,
					&n_bytes_in_packed_pixvals,
					&n_compressed_bytes_1,
					&n_compressed_bytes_2,
					run_metrics,
					copy_compressed_frame_to_return_buffer,    // copy_compressed_frame_to_return_buffer = 1, if 1 returns compressed_frame and compressed_frame_length, otherwise does not
					compressed_frame,
					&compressed_frame_length
				);

				if (available_buffer_space >= compressed_frame_length) {

					// copy frame to buffer
					for (i = 0; i < compressed_frame_length; i++) {
						pBuffer[RCT_Buffer_Fill_Position + i] = compressed_frame[i];
					}
					RCT_Buffer_Fill_Position += compressed_frame_length;
					available_buffer_space -= compressed_frame_length;
					//recode_print("RCT %d: Copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
					count++;

				} else {

					// offload data to disk
					//recode_print("RCT %d: buffer is full... clearing\n", process_id);
					fwrite(pBuffer, sizeof(char), RCT_Buffer_Fill_Position, partFile);
					fflush(partFile);
					nWrites++;
					available_buffer_space = BUFFER_LENGTH;
					RCT_Buffer_Fill_Position = 0;
					//recode_print("RCT %d: done.\n", process_id);

					// copy the last received frame that was never written
					for (i = 0; i < compressed_frame_length; i++) {
						pBuffer[RCT_Buffer_Fill_Position + i] = compressed_frame[i];
					}	
					RCT_Buffer_Fill_Position += compressed_frame_length;
					available_buffer_space -= compressed_frame_length;
					count++;
					//recode_print("RCT %d: Copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);

				}
				
			}

			count = 0;
			run_offset += n_frames_in_chunk;
			num_frames_in_part += num_frames_to_process;
			RCT_STATE_INDICATORS[process_id] = 0;

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
	//uint64_t BUFFER_LENGTH = params->num_frames*params->num_rows*params->num_cols;
	uint64_t BUFFER_LENGTH = 5*params->num_rows*params->num_cols;
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

void mergePartFiles(InitParams *init_params, InputParams *params, RCHeader *header) {

	const char *compressed_filename = makeCompressedFilename(init_params->image_filename, params->reduction_level);
	char** part_filenames = (char**)malloc((params->num_threads) * sizeof(char*));
	int part_num;
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		part_filenames[part_num] = makePartFilename(part_num, init_params->run_name, params->reduction_level);
	}
	merge_RC1_Parts(init_params->output_directory, part_filenames, header, params, compressed_filename);
}

int start_zmq_server(InitParams *init_params, InputParams *params, RCHeader *header) {

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
			//mergePartFiles(init_params, params, header);
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

void _simulate_zmq_run(InitParams *init_params, InputParams *params, RCHeader *header) {

	int iters = 3;
	for (int i = 0; i < iters; i++) {
		printf("ReCoDe Server: Received Request %d, processing: %d\n", i, i);
		reset_run_stats();
		//IMAGE_FILENAME = "R:\\003.seq";
		//IMAGE_FILENAME = "D:\\cbis\\GitHub\\ReCoDe\\scratch\\temp\\003.seq";
		IMAGE_FILENAME = "D:\\cbis\\images\\SequenceBlocks\\17-21-33.138_Dark_Ref_3.seq";
		set_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1);
		Sleep(3000);
		while (get_RCT_States(RCT_STATE_INDICATORS, params->num_threads, 1) != 0) {
			Sleep(100);          //  Wait for RCT Threads to finish
		}
		//mergePartFiles(init_params, params, header);
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
	/*====================Testing========================*/

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
	const char *image_name = getFilenameFromPath(init_params->image_filename);
	const char *dark_name = getFilenameFromPath(init_params->dark_filename);
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
		uint64_t sz_darkBuffer = (uint64_t)header->nx * (uint64_t)header->ny * (uint64_t)params->num_dark_frames * (uint64_t)sizeof(uint16_t);
		printf("Size of uint64_t = %d\n", sizeof(uint64_t));
		//uint64_t sz_darkBuffer = header->nx * header->ny * params->num_dark_frames * 2;
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
			//start_zmq_server (init_params, params, header);
			_simulate_zmq_run(init_params, params, header);
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

