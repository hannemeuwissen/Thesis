/**
 * @file decomp1d.c
 * @author Hanne Meuwissen (22307813)
 * @brief Code file for thesis at Trinity College Dublin.
 * @version 0.1
 * @date 2023-05-27
 */
#include <stdlib.h>
#include <stdio.h>

/**
 * @brief Function that decomposes an array of doubles of size n across processors.
 * @param[in] n The length of the array.
 * @param[in] p The number of processors.
 * @param[in] myid The rank of the processor.
 * @param[in] s Pointer to store the address to the startpoint of the rank in the array.
 * @param[in] e Pointer to store the address to the endpoint of the rank in the array.
 * @return The function returns 0 on success, -1 otherwise.
 */
int decomp1d(int n, int p, int myid, int *s, int *e){
    if((n<1) || (p<1) || (p>n)){
        perror("Invalid input.");
        return(-1);
    }
    int n_elements = n/p;
    int remainder = n%p;
    *e = 0;
    for(int i=0;i<=myid;i++){
        *s = *e + 1;
        *e += ((i<remainder) ? (n_elements + 1) : (n_elements)); 
    }
    *s-=1;
    *e-=1;
    return 0;
}