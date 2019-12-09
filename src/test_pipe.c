#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>


#include "tlrng.h"

#define BUFFSIZE 10

int main ()
{

        fprintf(stdout,"writing gzx file\n");
        char* buf = malloc(BUFFSIZE);
        char* wb;
        int i;
        gzFile fp;
        fp = gzopen("testzlib.gz", "wb");



        int index = 0;
        int c;

        wb = buf;
        for(i = 41; i < 100; i++){
                fprintf(stdout,"Space: %d\n", BUFFSIZE-index);
                c = snprintf(buf+index,BUFFSIZE-index, "%d,",i);

                if(c > BUFFSIZE-index){
                        fprintf(stdout,"No space!\n");
                        fprintf(stdout,"%s\n",buf);
                        gzwrite(fp, wb, index);
                        index = 0;
                }else{
                        index = index+c;
                }
        }
        if(index){
                gzwrite(fp, wb, index);
        }
        //gzwrite(fp, buf, len);
        //gzprintf(fp, "Hello, %s!\n", "world");
        gzclose(fp);

        free(buf);
        return EXIT_SUCCESS;
}
