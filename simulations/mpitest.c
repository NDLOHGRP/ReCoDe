/**************************************************************************************************************
to compile:    mpicc mpitest.c -L../zlib-1.2.11 -lz -L../bzip2-1.0.6 -lbz2 -L../lz4-1.8.2/lib -llz4 -O3 -o mpitest -lm -fopenmp
to run:        mpirun -n 2 -ppn 1 ./mpitest
               

-lz:     links zlib
-lbz2:   links bzip2
-lm:     links math.h NOTE: it is important to keep -lm as the last option
-o:      write build output to file
-O3:     optimizer flag
**************************************************************************************************************/



#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "mpi.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "zlib.h"
#include "bzlib.h"
//#include "lzma.h"
//#include "snappy-c.h"
#include "lz4.h"

#include <dirent.h>

#include "recodefs.h"
#include "recode_utils.h"
#include "commons.h"
#include "compressor.h"
#include "fileutils.h"
#include "logger.h"
#include "mrc_header.h"
#include "argparser.h"
#include "recode_header.h"

#include "L1_next.h"

#define    MASTER		                0
#define    RC_WRITE_MODE_NO_RETURN      0        // reduceCompressFrame_L1 does not write to return buffer
#define    RC_WRITE_MODE_RETURN_ONLY    1        // reduceCompressFrame_L1 writes to return buffer but the data is not copied to OLT buffer
#define    RC_WRITE_MODE_RETURN_TO_OLT  2        // data is copied to OLT buffer but not written to file
#define    RC_WRITE_MODE_FILE_WRITE     3        // data is written to file


// Global Variables
unsigned int    WAIT_MILLISECONDS    = 1;
unsigned int    ALL_RCT_DONE         = 0;
unsigned long   BUFFER_LENGTH;

// reduce-compress thread
unsigned long runRCT_L1 (   uint8_t rc_operation_mode,
                            uint8_t rc_write_mode,
                            unsigned short* frameBuffer, 
                            unsigned short* darkFrame, 
                            DataSize h, 
                            int process_id, 
                            int MPI_PID_Offset,
                            int num_frames_to_process, 
                            unsigned char *RCT_FULL_IND, 
                            unsigned long *RCT_Buffer_Fill_Position, 
                            unsigned char *pBuffer,
                            float *compression_times) {

    unsigned int i;

    unsigned long count = 0;

    unsigned long wait_time = 0;

    unsigned long available_buffer_space = BUFFER_LENGTH;

    unsigned long frame_offset = process_id * num_frames_to_process;
    
    
    uint8_t copy_compressed_frame_to_return_buffer;
    if (rc_write_mode > RC_WRITE_MODE_NO_RETURN) {
        copy_compressed_frame_to_return_buffer = 1;
    } else {
        copy_compressed_frame_to_return_buffer = 0;
    }


    uint32_t n_bytes_in_packed_pixvals, n_compressed_bytes_1, n_compressed_bytes_2, compressed_frame_length;
        
    // create data buffers
    uint32_t n_pixels_in_frame             = h.nx * h.ny;
    uint32_t n_bytes_in_binary_image   = ceil(n_pixels_in_frame / 8.0);
    uint16_t *pixvals                     	  = (uint16_t*) malloc (n_pixels_in_frame * sizeof(uint16_t));            // allocated max space - where all pixels are fg pixels, must calloc to init to 0
    uint8_t  *binaryImage                 = (uint8_t*) malloc (n_bytes_in_binary_image * sizeof(uint8_t));    // every bit in binaryImage has to be reset for every frame in the loop below
    uint8_t  *packedPixvals             = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));    
    uint8_t  *compressedPixvals         = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
    uint8_t  *compressedBinaryImage     = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
    uint8_t  *compressedPackedPixvals    = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
    uint8_t  *compressed_frame            = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
    
    int pid = MPI_PID_Offset + process_id;
    const char* filename = makePartFilename(pid, "test_11", 1);
	FILE *partFile = fopen(concat("out/", filename), "wb");
	//partFile = fopen(concat("../../../../../../scratch/loh/abhik/code_data/DataCompression/on_the_fly_output/",filename), "wb");
    //partFile = fopen(concat("/gpfs0/scratch/loh/abhik/code_data/DataCompression/on_the_fly_output/", filename), "wb");
    partFile = fopen(concat("/mnt/cbis/scratch-test/loh/abhik/code_data/DataCompression/on_the_fly_output/", filename), "wb");
	int FD = fileno(partFile);
	fwrite (&pid, sizeof(char), 1, partFile);        // write process_id to partfile
    
    unsigned long nWrites = 0;
    while(1) {

        if (count == num_frames_to_process) {
            *RCT_FULL_IND = 3;                            // this RCT is done, indicate OLT to close RCT's file
            
            fclose(partFile);
            
            free(pixvals);
            free(binaryImage);
            free(packedPixvals);
            free(compressedPixvals);
            free(compressedBinaryImage);
            free(compressedPackedPixvals);
            free(compressed_frame);
            recode_print("RCT %d: Done (total wait time = %u milliseconds).\n", process_id, wait_time * WAIT_MILLISECONDS);
            //return wait_time;
            return nWrites;
        }

        if (*RCT_FULL_IND == 2) {                        // buffer has been cleared by OLT, reset indicators
        
            available_buffer_space = BUFFER_LENGTH;
            *RCT_Buffer_Fill_Position = 0;
            *RCT_FULL_IND = 0;
            recode_print("RCT %d: OLT cleared buffer.\n", process_id);
            
            // copy the last received frame that was never written
            if (rc_write_mode > RC_WRITE_MODE_RETURN_ONLY) {
                for (i = 0; i<compressed_frame_length; i++) {
                    pBuffer[(*RCT_Buffer_Fill_Position) + i] = compressed_frame[i];
                }
                //recode_print("Copying to OLT Buffer.\n");
                
                recode_print("RCT %d \t Waited Milliseconds: \t%d.\n", process_id, WAIT_MILLISECONDS*wait_time);
                wait_time = 0;
            }
            
            (*RCT_Buffer_Fill_Position) = (*RCT_Buffer_Fill_Position) + compressed_frame_length;
            available_buffer_space -= compressed_frame_length;
            recode_print("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
            
        }

        if (*RCT_FULL_IND == 0) {                        // buffer has space, continue writing
            
            reduceCompressFrame_L1 (process_id,
                                    frameBuffer, 
                                    darkFrame, 
                                    frame_offset + count,                // absolute frame index
                                    frame_offset + count,                 // frame index relative to start of this thread
                                    h, 
                                    rc_operation_mode,                    // rc_operation_mode: 1 = reduce and compress
                                    0,                                    // epsilon_s
                                    12,                                    // bit_depth
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
                                    compression_times,
                                    copy_compressed_frame_to_return_buffer,    // copy_compressed_frame_to_return_buffer
                                    compressed_frame,
                                    &compressed_frame_length);
                                            
            //recode_print("RCT %d: writing frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);

            // check if buffer has enough space to write the next set of compressed bytes
            if (available_buffer_space >= compressed_frame_length) {
                
                // write to designated buffer
                // memcpy (&pBuffer[RCT_Buffer_Fill_Position], compressed_frame, compressed_frame_length);
                if (rc_write_mode > RC_WRITE_MODE_RETURN_ONLY) {
                    
                    for (i = 0; i<compressed_frame_length; i++) {
                        pBuffer[(*RCT_Buffer_Fill_Position) + i] = compressed_frame[i];
                    }
                    //recode_print("Copying to OLT Buffer.\n");
                }
                
                (*RCT_Buffer_Fill_Position) = (*RCT_Buffer_Fill_Position) + compressed_frame_length;
                available_buffer_space -= compressed_frame_length;
                recode_print("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
                count++;
                    
            } else {

                *RCT_FULL_IND = 1;            // buffer is full, flip indicator and wait for OLT to clear
                recode_print("RCT %d: buffer is full, waiting for OLT to clear.\n", process_id);

            }

        } 

        if (*RCT_FULL_IND == 1) {
            /*
            fwrite (pBuffer, sizeof(char), (*RCT_Buffer_Fill_Position), partFile); 
            fdatasync(FD);
            nWrites++;
            */
            *RCT_FULL_IND = 2;
            
            //wait_time++;
            //usleep(WAIT_MILLISECONDS*1000);
        }
        

    }

}


void recode_on_the_fly (int level, 
                        uint8_t rc_operation_mode,        // rc_operation_mode: 1 = reduce and compress
                        uint8_t rc_write_mode,
                        unsigned short* frameBuffer, 
                        unsigned short* darkFrame, 
                        int NT_RCT,
						int NT_OLT,
                        unsigned long BUFFER_LENGTH_IN, 
                        int nImageFrames, 
                        int doMemTracking, 
                        float dose_rate,
                        int MPI_PID,
						int replicate) {

    int i;

    srand((unsigned)time(NULL));

    BUFFER_LENGTH = BUFFER_LENGTH_IN;

    

    /*
    ========================================================================================== 
    Create RCT full indicators
    ==========================================================================================
    */
    unsigned char *RCT_FULL_INDs = (unsigned char*)calloc(NT_RCT, sizeof(char));
    for (i = 0; i<NT_RCT; i++) {
        RCT_FULL_INDs[i] = 0;
    }


    /*
    ========================================================================================== 
    Create array for holding buffer fillup statuses
    ==========================================================================================
    */
    unsigned long *RCT_Buffer_Fill_Positions = (unsigned long*)calloc(NT_RCT, sizeof(long));
    for (i = 0; i<NT_RCT; i++) {
        RCT_Buffer_Fill_Positions[i] = 0;
    }

    /*
    ========================================================================================== 
    Create RCT buffers
    ==========================================================================================
    */
    unsigned char **pBuffers;
    pBuffers = (unsigned char**) malloc (NT_RCT * sizeof(unsigned char *));
    for (i = 0; i < NT_RCT; i++) {
          pBuffers[i] = (char *)malloc(BUFFER_LENGTH * sizeof(char)); // make actual arrays
    }
    
	/*
    ========================================================================================== 
    Create RCT-OLT association
    ==========================================================================================
    */
    uint8_t *Connected_OLTs = (uint8_t*)calloc(NT_RCT, sizeof(uint8_t));
	int n_RCTs_per_OLT = (int)ceil((NT_RCT*1.0)/(NT_OLT*1.0));
	int Connected_OLT_ID = -1;
    for (i = 0; i<NT_RCT; i++) {
		if (i%n_RCTs_per_OLT == 0) {
			Connected_OLT_ID++;
		}
        Connected_OLTs[i] = Connected_OLT_ID;
		recode_print("Connected-OLT: %d\n", Connected_OLT_ID);
    }

    /*
    ========================================================================================== 
    Do multithreaded reduction-compression
    ==========================================================================================
    */
    unsigned long nFramesPerThread = nImageFrames/NT_RCT;
    DataSize h = { 4096, 512, nFramesPerThread, 16 };

    unsigned long *num_writes = (unsigned long*)calloc(NT_RCT, sizeof(long));
    for (i = 0; i<NT_RCT; i++) {
        num_writes[i] = 0;
    }
    
    
    float *compress_times_sizes_best_speed = (float*) malloc (NT_RCT * 2 * sizeof(float));
    
    float total_time;
    double start = omp_get_wtime();
    
    if (doMemTracking == 0) {
    
        #pragma omp parallel for num_threads(NT_RCT + NT_OLT)
        for (i = 0; i < (NT_RCT + NT_OLT); i++) {

            if (i < NT_OLT) {
                
                //runOLT (NT_RCT, Connected_OLTs, i, rc_write_mode, &RCT_FULL_INDs, &RCT_Buffer_Fill_Positions, pBuffers);
                
            } else {

                if (level == 1) {

                    num_writes[i] = runRCT_L1 ( rc_operation_mode, 
                                    rc_write_mode,
                                    frameBuffer, 
                                    darkFrame, 
                                    h, 
                                    i-NT_OLT, 
                                    MPI_PID*NT_RCT,
                                    nFramesPerThread, 
                                    &RCT_FULL_INDs[i-NT_OLT], 
                                    &RCT_Buffer_Fill_Positions[i-NT_OLT], 
                                    pBuffers[i-NT_OLT], 
                                    compress_times_sizes_best_speed);

                }

            }

        }

        total_time = omp_get_wtime()-start;
		printf ("OLT Threads = %d, RCT Threads = %d\n", NT_OLT, NT_RCT); 
        printf ("Time:%f\n", total_time); 
        

    }


    /*
    ========================================================================================== 
    Log Results
    ========================================================================================== 
    */
    float avg_wait_time_per_thread = 0.0;
    float avg_write_events_per_thread = 0.0;
    for (i = 0; i<NT_RCT; i++) {
        avg_wait_time_per_thread += compress_times_sizes_best_speed[i];
        avg_write_events_per_thread += num_writes[i]*1.0;
    }
    avg_wait_time_per_thread /= NT_RCT;
    avg_wait_time_per_thread /= 1000;
    avg_write_events_per_thread /= (NT_RCT*1.0);
    printf ("Avg. Wait Time:\t%f (%f perc. of total time)\n", avg_wait_time_per_thread, 
                                                        avg_wait_time_per_thread/total_time);

    FILE *logfile2 = fopen("recode_on_the_fly_results", "a");
    
    fprintf(logfile2, "L%d \t %.1f \t %d \t %lu \t %lu \t %.4f \t %d \t %d \t %d \t %d \t %d \t %.3f\n", 
                        level, dose_rate, NT_RCT, nImageFrames, BUFFER_LENGTH_IN, total_time,
                        rc_operation_mode, rc_write_mode, NT_OLT, MPI_PID, replicate, avg_write_events_per_thread);
    fclose(logfile2);



    /*
    ========================================================================================== 
    Cleanup
    ==========================================================================================
    */
    free(RCT_FULL_INDs);
    free(RCT_Buffer_Fill_Positions);
    recode_print("2\n");
    for (i = 0; i < NT_RCT; i++) {
          free(pBuffers[i]);
    }
    //free(pBuffers);
	free(Connected_OLTs);
    free(compress_times_sizes_best_speed);
    free(num_writes);
    recode_print("3\n");
        
    return;

}

int main(int argc, char *argv[]) {
    
    int nNodes, taskid;
    MPI_Status status;

	/***** Initializations *****/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	printf("MPI task %d has started...\n", taskid);
 
    if (taskid == MASTER) {
        printf("Number of MPI Processes = %d\n", nNodes);
    }
    
    // Get the name of the processor
    char processor_name[255];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Task %d running on: %s\n", nNodes, processor_name);
    
    /*
    char *ptr;
    taskid = strtol(argv[1], &ptr, 10);
    */
    
    double dose_rates[4] = {1.6};
    
    int dr;
    for (dr = 0; dr < 4; dr++) {

        double dose_rate = dose_rates[dr];
        
        /*
        ========================================================================================== 
        Load Image Data
        ==========================================================================================
        */
        const char* imageFile;
        int nImageFrames;
        if (dose_rate == 0.8) {
            //imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-06-23.598.bin";
            //imageFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-06-23.598.bin";
            //imageFile = "/gpfs0/scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-06-23.598";
            imageFile = "/mnt/cbis/scratch-test/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-06-23.598";
            nImageFrames = 700*20;
        } else if (dose_rate == 1.6) {
            //imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-13-59.436.bin";
            //imageFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-13-59.436.bin";
            //imageFile = "/gpfs0/scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-13-59.436";
            imageFile = "/mnt/cbis/scratch-test/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-13-59.436";
            nImageFrames = 700*20;
        } else if (dose_rate == 3.2) {
            //imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-12-04.498.bin";
            //imageFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-12-04.498.bin";
            //imageFile = "/gpfs0/scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-12-04.498";
            imageFile = "/mnt/cbis/scratch-test/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-12-04.498";
            nImageFrames = 700*20;
        } else if (dose_rate == 6.4) {
            //imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-19-00.764.bin";
            //imageFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-19-00.764.bin";
            //imageFile = "/gpfs0/scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-19-00.764";
            imageFile = "/mnt/cbis/scratch-test/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-19-00.764";
            nImageFrames = 700*20;
        }
        
        char* taskid_str = malloc(sizeof(char)*1);
        sprintf(taskid_str, "_%d", taskid);
        DataSize b = { 4096, 512, nImageFrames, 16 };
        uint16_t *frameBuffer = (uint16_t *)malloc(b.nx * b.ny * b.nz * sizeof(uint16_t));
        
        if (taskid == MASTER) {
            const char* process_imageFile = concat(imageFile, ".bin");
            printf("%s\n", process_imageFile);
            loadData(process_imageFile, frameBuffer, 0, nImageFrames, b);
        } else {
            const char* process_imageFile = concat(concat(imageFile, taskid_str), ".bin");
            printf("%s\n", process_imageFile);
            loadData(process_imageFile, frameBuffer, 0, nImageFrames, b);
        }
        
        
        recode_print("Node %d: Image data loaded.\n", taskid);
        
        
        /*
        ========================================================================================== 
        Dark Noise Estimation
        ==========================================================================================
        */
        uint16_t *darkFrame = (uint16_t *)calloc(4096*512, sizeof(uint16_t));
        
        //char* darkNoiseFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/Dark_Frame_12-23-00.232.bin";
        //char* darkNoiseFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/Dark_Frame_12-23-00.232.bin";
        //char* darkNoiseFile = "/gpfs0/scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/Dark_Frame_12-23-00.232";
        char* darkNoiseFile = "/mnt/cbis/scratch-test/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/Dark_Frame_12-23-00.232";
        char* process_darkNoiseFile = concat(concat(darkNoiseFile, taskid_str), ".bin");
        DataSize d1 = { 4096, 512, 1, 12 };
        printf("%s\n", process_darkNoiseFile);
        loadData (process_darkNoiseFile, darkFrame, 0, 1, d1);
        recode_print("Node %d: Pre-computed dark noise estimates loaded.\n", taskid);
        
        
        if (taskid == MASTER) {
            
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
            fprintf(logfile, "Data NOT Saved (running on 2 Nodes): Kylin\n");
            fprintf(logfile, "=======================================\n");
            fprintf(logfile, "L \t DR \t NT \t n_Fr \t B_Sz \t T \t R/C \t WMODE \t OLT \t MPID \t NWE\n");
            fclose(logfile);
            
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        uint8_t reps;
        uint8_t level;
        uint8_t rc_operation_mode;            // rc_operation_mode: 1 = reduce and compress
        uint8_t rc_write_mode;
        uint8_t rct_threads;
		uint8_t olt_threads;
        
        for (level = 1; level < 2; level += 3) {
            for (rct_threads = 1; rct_threads < 29; rct_threads = rct_threads + 2) { 
				for (olt_threads = 0; olt_threads < 1; olt_threads++) {
					for (rc_operation_mode = 1; rc_operation_mode < 2; rc_operation_mode++) {
						for (rc_write_mode = 3; rc_write_mode < 4; rc_write_mode++) {
							for (reps = 0; reps < 5; reps++) { 
                            
								recode_on_the_fly ( level, 
													rc_operation_mode,
													rc_write_mode,
													frameBuffer, 
													darkFrame, 
													rct_threads, 
													olt_threads,
													1024*1024*10, 
													nImageFrames, 
													0, 
													dose_rate,
                                                    taskid,
													reps);
                                                    
                                                    
                                MPI_Barrier(MPI_COMM_WORLD);
                                
							}
						}
					}
				}
			}
        }

        
        
        recode_print("Node %d: Crossed Barrier.\n", taskid);
        
        /*
        ========================================================================================== 
        Cleanup
        ==========================================================================================
        */
        printf("0\n");
        free(darkFrame);
        printf("1\n");
        free(frameBuffer);
        printf("2\n");
        
        printf("Node %d: Done.\n", taskid);
    }
 
    MPI_Finalize();

}




