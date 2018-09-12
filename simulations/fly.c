/**************************************************************************************************************
to compile:    gcc fly.c -lz -lbz2 -llzma -lsnappy -llz4 -O3 -o fly -lm -fopenmp
to compile with profiling: kinst-ompp-papi gcc fly.c -lz -O3 -o fly -lm -fopenmp

to run:        ./fly

-lz:     links zlib
-lbz2:     links bzip2
-lm:     links math.h NOTE: it is important to keep -lm as the last option
-o:         write build output to file
-O3:     optimizer flag
**************************************************************************************************************/

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "zlib.h"
#include <bzlib.h>
#include <lzma.h>
#include "snappy-c.h"
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
#include "L4.h"




#define    RC_WRITE_MODE_NO_RETURN                0        // reduceCompressFrame_L1 does not write to return buffer
#define    RC_WRITE_MODE_RETURN_ONLY            1        // reduceCompressFrame_L1 writes to return buffer but the data is not copied to OLT buffer
#define    RC_WRITE_MODE_RETURN_TO_OLT            2        // data is copied to OLT buffer but not written to file
#define    RC_WRITE_MODE_FILE_WRITE            3        // data is written to file


// Global Variables
unsigned int     WAIT_MILLISECONDS    = 1;
unsigned int    ALL_RCT_DONE        = 0;
unsigned long     BUFFER_LENGTH;


//unsigned short* frameBuffer;
//unsigned short* darkFrame;


//stub to simulate reducer-compressor
int do_L1_Reduce_Compress_In_Memory (unsigned char** compressedData, unsigned long *nCompressedSize) {

    clock_t p_start = clock();

    int r = rand() % 1000;
    *nCompressedSize = 220201 + r;

    *compressedData = (unsigned char*)malloc((*nCompressedSize)*sizeof(char));
    //memcpy (*compressedData, frameBuffer, *nCompressedSize);
    unsigned long i;
    unsigned char count = 0;
    for (i = 0; i<(*nCompressedSize); i++) {
        (*compressedData)[i] = count++;
        if (count == 255) {
            count = 0;
        }
    }

    clock_t p_end = clock();
    float process_time = (p_end - p_start) * 1000.0 / CLOCKS_PER_SEC;

    //printf("process_time = %f\n", process_time);
    usleep(110*1000);

    return;
}


// memory-tracking thread
void runMTT (int NT, unsigned long **RCT_Buffer_Fill_Positions, unsigned long MAX_POLLING_ITERS, int POLLING_TIME_STEP) {
    
    unsigned long *maps = (unsigned long*)calloc(MAX_POLLING_ITERS, sizeof(long));
    
    unsigned long i;
    for (i = 0; i<MAX_POLLING_ITERS; i++) {
        maps[i] = 0;
    }
    
    printf("MMT: Tracking Started.\n");
    unsigned long count = 0;
    while (count < MAX_POLLING_ITERS && ALL_RCT_DONE == 0) {
        
        maps[count] = 0;
        for (i = 0; i<NT; i++) {
            maps[count] += (*RCT_Buffer_Fill_Positions)[i];
        }
        
        usleep(POLLING_TIME_STEP*1000);
        count++;
    }
    
    FILE *fp = fopen("Memmap.txt", "w");
    for (i = 0; i<count; i++) {
        fprintf (fp, "%d\t%u\n", i*POLLING_TIME_STEP, maps[i]);
    }
    fclose(fp);
    
    printf("MMT: Tracking Stopped.\n");
}



// off-loader thread
void runOLT (int NT, 
             uint8_t rc_write_mode,
             unsigned char **RCT_FULL_INDs, 
             unsigned long **RCT_Buffer_Fill_Positions, 
             unsigned char **pBuffers) {

    int i;

    FILE *partFiles[NT-1];
    int FDs[NT-1];
    for (i = 0; i < NT; i++) {
        const char* filename = makePartFilename(i, "test_11", 1);
        //partFiles[i] = fopen(concat("out/",filename), "wb");
        partFiles[i] = fopen(concat("../../../../../../scratch/loh/abhik/code_data/DataCompression/on_the_fly_output/",filename), "wb");
        FDs[i] = fileno(partFiles[i]);
        fwrite (&i, sizeof(char), 1, partFiles[i]);        // write process_id to partfile
    }
    ///printf("OLT Ready.\n");

    int num_open_threads = NT;
    while (1) {

        for (i = 0; i < NT; i++) {
            
            if ((*RCT_FULL_INDs)[i] == 1) {                // RCT indicates buffer is full, off load data to file
                ///printf("OLT: RCT %d's buffer is full, trying to clear %u bytes to disk.\n", i, (*RCT_Buffer_Fill_Positions)[i]);
                
                if (rc_write_mode > RC_WRITE_MODE_RETURN_TO_OLT) {
                    fwrite (pBuffers[i], sizeof(char), 
                            (*RCT_Buffer_Fill_Positions)[i], 
                            partFiles[i]);                        // off-load pBuffers[i] to disk
                    fdatasync(FDs[i]);
                    ///printf("Writing to disk.\n");
                }
                
                (*RCT_FULL_INDs)[i] = 2;                    // indicate RCT buffer has been off-loaded
                ///printf("OLT: RCT %d is cleared as requested.\n", i);
            }

            if ((*RCT_FULL_INDs)[i] == 3) {                // RCT indicates it is done, close RCT's file
                
                
                ///printf("OLT: RCT %d indicates it is done. Flushing available buffer data to file.\n", i);
                
                if (rc_write_mode > RC_WRITE_MODE_RETURN_TO_OLT) {
                    fwrite (pBuffers[i], sizeof(char), 
                            (*RCT_Buffer_Fill_Positions)[i], 
                            partFiles[i]);                // off-load pBuffers[i] to disk    
                    fdatasync(FDs[i]);
                    ///printf("Writing to disk.\n");
                }
                
                //fclose(partFiles[i]);
                num_open_threads--;
                (*RCT_FULL_INDs)[i] = 4;                // RCT is now closed
                ///printf("OLT: RCT %d is closed.\n", i);
            }

        }

        if (num_open_threads == 0) {                    // All RCTs are done, exit
        
            for (i = 0; i < NT; i++) {
                fclose(partFiles[i]);
            }
        
            ALL_RCT_DONE = 1;
            ///printf("OLT: All RCTs are done.\n");
            return;
        }

    }

}

// reduce-compress thread
unsigned long runRCT_L1 (    uint8_t rc_operation_mode,
                            uint8_t rc_write_mode,
                            unsigned short* frameBuffer, 
                            unsigned short* darkFrame, 
                            DataSize h, 
                            int process_id, 
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
    
    
    while(1) {

        if (count == num_frames_to_process) {
            *RCT_FULL_IND = 3;                            // this RCT is done, indicate OLT to close RCT's file
            free(pixvals);
            free(binaryImage);
            free(packedPixvals);
            free(compressedPixvals);
            free(compressedBinaryImage);
            free(compressedPackedPixvals);
            free(compressed_frame);
            ///printf("RCT %d: Done (total wait time = %u milliseconds).\n", process_id, wait_time * WAIT_MILLISECONDS);
            return wait_time;
        }

        if (*RCT_FULL_IND == 2) {                        // buffer has been cleared by OLT, reset indicators
        
            available_buffer_space = BUFFER_LENGTH;
            *RCT_Buffer_Fill_Position = 0;
            *RCT_FULL_IND = 0;
            //printf("RCT %d: OLT cleared buffer.\n", process_id);
            
            // copy the last received frame that was never written
            if (rc_write_mode > RC_WRITE_MODE_RETURN_ONLY) {
                for (i = 0; i<compressed_frame_length; i++) {
                    pBuffer[(*RCT_Buffer_Fill_Position) + i] = compressed_frame[i];
                }
                ///printf("Copying to OLT Buffer.\n");
                
                ///printf("RCT %d \t Waited Milliseconds: \t%d.\n", process_id, WAIT_MILLISECONDS*wait_time);
                wait_time = 0;
            }
            
            (*RCT_Buffer_Fill_Position) = (*RCT_Buffer_Fill_Position) + compressed_frame_length;
            available_buffer_space -= compressed_frame_length;
            //printf("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
            
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
                                    1,                                    // compression_scheme: 1 = gzip
                                    0,                                    // compression_level: 0 = fastest
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
                                            
                                            
            //printf("RCT %d: writing frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);

            // check if buffer has enough space to write the next set of compressed bytes
            if (available_buffer_space >= compressed_frame_length) {
                
                // write to designated buffer
                // memcpy (&pBuffer[RCT_Buffer_Fill_Position], compressed_frame, compressed_frame_length);
                if (rc_write_mode > RC_WRITE_MODE_RETURN_ONLY) {
                    for (i = 0; i<compressed_frame_length; i++) {
                        pBuffer[(*RCT_Buffer_Fill_Position) + i] = compressed_frame[i];
                    }
                    ///printf("Copying to OLT Buffer.\n");
                }
                
                (*RCT_Buffer_Fill_Position) = (*RCT_Buffer_Fill_Position) + compressed_frame_length;
                available_buffer_space -= compressed_frame_length;
                //printf("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
                count++;
                    
            } else {

                *RCT_FULL_IND = 1;            // buffer is full, flip indicator and wait for OLT to clear
                ///printf("RCT %d: buffer is full, waiting for OLT to clear.\n", process_id);

            }

        } 

        if (*RCT_FULL_IND == 1) {
            wait_time++;
            usleep(WAIT_MILLISECONDS*1000);
        }
        

    }

}


// reduce-compress thread
unsigned long runRCT_L4 (    uint8_t rc_operation_mode,
                            uint8_t rc_write_mode,
                            unsigned short* frameBuffer, 
                            unsigned short* darkFrame, 
                            DataSize h, 
                            int process_id, 
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
    

    uint32_t n_compressed_bytes, compressed_frame_length;
    
    
    uint8_t copy_compressed_frame_to_return_buffer;
    if (rc_write_mode > RC_WRITE_MODE_NO_RETURN) {
        copy_compressed_frame_to_return_buffer = 1;
    } else {
        copy_compressed_frame_to_return_buffer = 0;
    }
    
    uint32_t n_pixels_in_frame             = h.nx * h.ny;
    uint32_t n_bytes_in_binary_image     = ceil(n_pixels_in_frame / 8.0);
    float      *x_coor                      = (float *) malloc  (n_pixels_in_frame * sizeof(float));
    float      *y_coor                      = (float *) malloc  (n_pixels_in_frame * sizeof(float));
    uint16_t *foregroundImage              = (uint16_t *)calloc(n_pixels_in_frame, sizeof(uint16_t));
    uint8_t  *centroidImage             = (uint8_t *)calloc (n_bytes_in_binary_image, 1);            // allocated max space - where all pixels are fg pixels, must calloc to init to 0
    uint8_t  *foregroundTernaryMap         = (uint8_t *)calloc (n_pixels_in_frame, sizeof(uint8_t));
    uint8_t  *compressedCentroidImage    = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
    uint8_t  *compressed_frame            = (uint8_t*) malloc (n_pixels_in_frame * sizeof(uint8_t));
            
    
    while(1) {

        if (count == num_frames_to_process) {
            *RCT_FULL_IND = 3;                            // this RCT is done, indicate OLT to close RCT's file
            free(x_coor);
            free(y_coor);
            free(foregroundImage);
            free(centroidImage);
            free(foregroundTernaryMap);
            free(compressedCentroidImage);
            free(compressed_frame);
            //printf("RCT %d: Done (total wait time = %u milliseconds).\n", process_id, wait_time * WAIT_MILLISECONDS);
            return wait_time;
        }

        if (*RCT_FULL_IND == 2) {                        // buffer has been cleared by OLT, reset indicators
        
            available_buffer_space = BUFFER_LENGTH;
            *RCT_Buffer_Fill_Position = 0;
            *RCT_FULL_IND = 0;
            //printf("RCT %d: OLT cleared buffer.\n", process_id);
            
            // copy the last received frame that was never written
            for (i = 0; i<compressed_frame_length; i++) {
                pBuffer[(*RCT_Buffer_Fill_Position) + i] = compressed_frame[i];
            }
                
            (*RCT_Buffer_Fill_Position) = (*RCT_Buffer_Fill_Position) + compressed_frame_length;
            available_buffer_space -= compressed_frame_length;
            //printf("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
            
        }

        if (*RCT_FULL_IND == 0) {                        // buffer has space, continue writing
            
            
            reduceCompressFrame_L4 (process_id,
                                    frameBuffer, 
                                    darkFrame, 
                                    frame_offset + count,                // absolute frame index
                                    frame_offset + count,                 // frame index relative to start of this thread
                                    h, 
                                    rc_operation_mode,                    // rc_operation_mode: 1 = reduce and compress
                                    0,                                    // epsilon_s
                                    1,                                    // compression_scheme: 1 = gzip
                                    0,                                    // compression_level: 0 = fastest
                                    x_coor, 
                                    y_coor, 
                                    foregroundImage, 
                                    centroidImage,
                                    foregroundTernaryMap,
                                    compressedCentroidImage,
                                    &n_compressed_bytes, 
                                    compression_times,
                                    copy_compressed_frame_to_return_buffer,    // copy_compressed_frame_to_return_buffer
                                    compressed_frame,
                                    &compressed_frame_length);
                                            
                                            
            //printf("RCT %d: writing frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);

            // check if buffer has enough space to write the next set of compressed bytes
            if (available_buffer_space >= compressed_frame_length) {
                
                // write to designated buffer
                // memcpy (&pBuffer[RCT_Buffer_Fill_Position], compressed_frame, compressed_frame_length);
                for (i = 0; i<compressed_frame_length; i++) {
                    pBuffer[(*RCT_Buffer_Fill_Position) + i] = compressed_frame[i];
                }
                
                (*RCT_Buffer_Fill_Position) = (*RCT_Buffer_Fill_Position) + compressed_frame_length;
                available_buffer_space -= compressed_frame_length;
                //printf("RCT %d: copied frame %d (%u bytes) to buffer.\n", process_id, frame_offset + count, compressed_frame_length);
                count++;
                    
            } else {

                *RCT_FULL_IND = 1;            // buffer is full, flip indicator and wait for OLT to clear
                //printf("RCT %d: buffer is full, waiting for OLT to clear.\n", process_id);

            }

        } 

        if (*RCT_FULL_IND == 1) {
            wait_time++;
            usleep(WAIT_MILLISECONDS*1000);
        }
        

    }

}



void recode_on_the_fly (int level, 
                        uint8_t rc_operation_mode,        // rc_operation_mode: 1 = reduce and compress
                        uint8_t rc_write_mode,
                        unsigned short* frameBuffer, 
                        unsigned short* darkFrame, 
                        int NT, 
                        unsigned long BUFFER_LENGTH_IN, 
                        int nImageFrames, 
                        int doMemTracking, 
                        float dose_rate) {

    int i;

    srand((unsigned)time(NULL));

    BUFFER_LENGTH = BUFFER_LENGTH_IN;

    

    /*
    ========================================================================================== 
    Create RCT full indicators
    ==========================================================================================
    */
    unsigned char *RCT_FULL_INDs = (unsigned char*)calloc(NT, sizeof(char));
    for (i = 0; i<NT; i++) {
        RCT_FULL_INDs[i] = 0;
    }


    /*
    ========================================================================================== 
    Create array for holding buffer fillup statuses
    ==========================================================================================
    */
    unsigned long *RCT_Buffer_Fill_Positions = (unsigned long*)calloc(NT, sizeof(long));
    for (i = 0; i<NT; i++) {
        RCT_Buffer_Fill_Positions[i] = 0;
    }

    /*
    ========================================================================================== 
    Create RCT buffers
    ==========================================================================================
    */
    unsigned char **pBuffers;
    pBuffers = (unsigned char**) malloc (NT * sizeof(unsigned char *));
    for (i = 0; i < NT; i++) {
          pBuffers[i] = (char *)malloc(BUFFER_LENGTH * sizeof(char)); // make actual arrays
    }
    

    /*
    ========================================================================================== 
    Do multithreaded reduction-compression
    ==========================================================================================
    */
    unsigned long nFramesPerThread = nImageFrames/NT;
    DataSize h = { 4096, 512, nFramesPerThread, 16 };

    /*
    unsigned long *wait_times = (unsigned long*)calloc(NT, sizeof(char));
    for (i = 0; i<NT; i++) {
        wait_times[i] = 0;
    }
    */
    
    float *compress_times_sizes_best_speed = (float*) malloc (NT * 2 * sizeof(float));
    
    float total_time;
    double start = omp_get_wtime();
    
    if (doMemTracking == 0) {
    
        #pragma omp parallel for num_threads(NT+1)
        for (i = 0; i <= NT; i++) {

            if (i == 0) {
                
                runOLT (NT, rc_write_mode, &RCT_FULL_INDs, &RCT_Buffer_Fill_Positions, pBuffers);
                
            } else {

                if (level == 1) {

                    runRCT_L1 (    rc_operation_mode, 
                                rc_write_mode,
                                frameBuffer, 
                                darkFrame, 
                                h, 
                                i-1, 
                                nFramesPerThread, 
                                &RCT_FULL_INDs[i-1], 
                                &RCT_Buffer_Fill_Positions[i-1], 
                                pBuffers[i-1], 
                                compress_times_sizes_best_speed);

                } else if (level == 4) {

                    runRCT_L4 (    rc_operation_mode, 
                                rc_write_mode,
                                frameBuffer, 
                                darkFrame, 
                                h, 
                                i-1, 
                                nFramesPerThread, 
                                &RCT_FULL_INDs[i-1], 
                                &RCT_Buffer_Fill_Positions[i-1], 
                                pBuffers[i-1], 
                                compress_times_sizes_best_speed);

                }

            }

        }

        total_time = omp_get_wtime()-start;
        printf("Time:\t\t%f \n", total_time); 
        

    } else {
        printf("1\n");
        
        #pragma omp parallel for num_threads(NT+2)
        for (i = 0; i <= NT+1; i++) {

            if (i == 0) {
                
                runOLT (NT, rc_write_mode, &RCT_FULL_INDs, &RCT_Buffer_Fill_Positions, pBuffers);
                
            } else if (i == 1) {

                runMTT (NT, &RCT_Buffer_Fill_Positions, 80, 250);
            
            } else {
                
                runRCT_L1 ( rc_operation_mode, 
                            rc_write_mode, 
                            frameBuffer, darkFrame, h, i-2, nFramesPerThread, 
                            &RCT_FULL_INDs[i-2], &RCT_Buffer_Fill_Positions[i-2], pBuffers[i-2], compress_times_sizes_best_speed);
                        
            }

        }

        total_time = omp_get_wtime()-start;
        printf("Threads = %d\tTime:\t\t%f \n", NT, total_time); 

        
    }


    /*
    ========================================================================================== 
    Log Results
    ========================================================================================== 
    */
    float avg_wait_time_per_thread = 0.0;
    for (i = 0; i<NT; i++) {
        avg_wait_time_per_thread += compress_times_sizes_best_speed[i];
    }
    avg_wait_time_per_thread /= NT;
    avg_wait_time_per_thread /= 1000;
    printf("Avg. Wait Time:\t%f (%f perc. of total time)\n", avg_wait_time_per_thread, 
                                                        avg_wait_time_per_thread/total_time);

    FILE *logfile2 = fopen("recode_on_the_fly_results", "a");
    
    fprintf(logfile2, "L%d \t %.1f \t %d \t %lu \t %lu \t %.4f \t %d \t %d\n", 
                        level, dose_rate, NT, nImageFrames, BUFFER_LENGTH_IN, total_time,
                        rc_operation_mode, rc_write_mode);
    fclose(logfile2);



    /*
    ========================================================================================== 
    Cleanup
    ==========================================================================================
    */
    free(RCT_FULL_INDs);
    free(RCT_Buffer_Fill_Positions);
    printf("2\n");
    for (i = 0; i < NT; i++) {
          free(pBuffers[i]);
    }
    //free(pBuffers);
    free(compress_times_sizes_best_speed);
    printf("3\n");
        
    return;

}


int main(int argc, char *argv[]) {


    int nImageFrames;
    int replicates = 1;
    
    
    double dose_rates[4] = {0.8, 1.6, 3.2, 6.4};
    
    int dr;
    for (dr = 0; dr < 4; dr++) {

        double dose_rate = dose_rates[dr];

        /*
        ========================================================================================== 
        Load Image Data
        ==========================================================================================
        */
        const char* imageFile;
        if (dose_rate == 0.8) {
            imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-06-23.598.bin";
            nImageFrames = 700*1;
        } else if (dose_rate == 1.6) {
            imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-13-59.436.bin";
            nImageFrames = 700*50;
        } else if (dose_rate == 3.2) {
            //imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-12-04.498.bin";
            imageFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-12-04.498.bin";
            nImageFrames = 700*20;
        } else if (dose_rate == 6.4) {
            imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-19-00.764.bin";
            nImageFrames = 700*20;
        }
        
        DataSize b = { 4096, 512, nImageFrames, 16 };
        uint16_t *frameBuffer = (uint16_t *)malloc(b.nx * b.ny * b.nz * sizeof(uint16_t));
        loadData(imageFile, frameBuffer, 0, nImageFrames, b);
        printf("Image data loaded.\n");
        
        /*
        ========================================================================================== 
        Replicate Image Data
        ==========================================================================================
        
        unsigned long r, f, row, col, count;
        unsigned long frame_sz = b.nx * b.ny;
        unsigned short* frameBuffer = (unsigned short*)malloc(b.nx * b.ny * b.nz * replicates * sizeof(short));
        
        count = 0;
        for (r = 0; r < replicates; r++) {
            for (f = 0; f < b.nz; f++) {
                for (row = 0; row < b.ny; row++) {
                    for (col = 0; col < b.nx; col++) {
                        frameBuffer[count++] = imageBuffer[f * frame_sz + row * b.nx + col];
                    }
                }
            }
        }
        
        nImageFrames = nImageFrames * replicates;
        b.nz = nImageFrames;
        printf("Image Data Replicated %d times.\n", replicates);
        */
        
        /*
        ========================================================================================== 
        Dark Noise Estimation
        ==========================================================================================
        */
        
        uint16_t *darkFrame = (uint16_t *)calloc(4096*512, sizeof(uint16_t));
        
        char* darkNoiseFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/Dark_Frame_12-23-00.232.bin";
        //char* darkNoiseFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/Dark_Frame_12-23-00.232.bin";
        DataSize d1 = { 4096, 512, 1, 12 };
        loadData (darkNoiseFile, darkFrame, 0, 1, d1);
        printf("Pre-computed dark noise estimates loaded.\n");
        
        
        /*
        int nDarkFrames = 10000;
        const char* darkFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-23-00.232.bin";
        DataSize d = { 4096, 512, nDarkFrames, 16 };
        uint16_t *darkBuffer = (uint16_t*)malloc(d.nx * d.ny * d.nz * sizeof(uint16_t));
        loadData(darkFile, darkBuffer, 0, nDarkFrames, d);
        printf("Dark data loaded.\n");
        getDarkMax (darkBuffer, d, darkFrame, 25);
        free(darkBuffer);
        printf("Dark noise estimation complete.\n");
        */
        
        /*
        int threads, mem;
        for (threads = 5; threads <=20; threads = threads + 5) {
            for (mem = 10; mem <= 50; mem = mem + 10) {
                printf("Threads = %d, Memory = %d MB\n", threads, mem);
                recode_on_the_fly (frameBuffer, darkFrame, threads, 1024*1024*mem, nImageFrames, 0);
            }
        }
        */

        
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
        fprintf(logfile, "Data Saved to Network Drive over 10GigE\n");
        fprintf(logfile, "=======================================\n");
        fprintf(logfile, "L \t DR \t NT \t n_Fr \t B_Sz \t T \t R/C \t WMODE\n");
        fclose(logfile);
        
        
        int reps;
        
        uint8_t level = 1;
        uint8_t rc_operation_mode = 0;            // rc_operation_mode: 1 = reduce and compress
        uint8_t rc_write_mode = 0;
        uint8_t threads = 0;
        
        for (level = 1; level < 2; level += 3) {
            for (threads = 5; threads < 51; threads = threads + 5) { 
                for (rc_operation_mode = 1; rc_operation_mode < 2; rc_operation_mode++) {
                    for (rc_write_mode = 3; rc_write_mode < 4; rc_write_mode++) {
                        for (reps = 0; reps < 5; reps++) { 
                            recode_on_the_fly (    level, 
                                                rc_operation_mode,
                                                rc_write_mode,
                                                frameBuffer, 
                                                darkFrame, 
                                                threads, 
                                                1024*1024*10, 
                                                nImageFrames, 
                                                0, 
                                                dose_rate);
                        }
                    }
                }
            }
        }
        
        
        /*
        for (threads = 25; threads < 26; threads = threads + 5) { 
            for (reps = 0; reps < 5; reps++) { 
                recode_on_the_fly (4, frameBuffer, darkFrame, threads, 1024*1024*10, nImageFrames, 0, dose_rate);
            }
        }
        */

        /*
        ========================================================================================== 
        Cleanup
        ==========================================================================================
        */
        printf("0\n");
        free(darkFrame);
        free(frameBuffer);
        //free(imageBuffer);
        printf("1\n");
        
    }
}
