

#include "io_hdf5.h"


#include "hdf5.h"
#include "hdf5_hl.h"

#include "tldevel.h"
#include "tlrng.h"


struct extract_plan{
        int** data;

};

int main(int argc, char *argv[])
{

        struct extract_plan* ep = NULL;
        hsize_t dims[2];
        int i,j;
        int file_id;
        MMALLOC(ep, sizeof(struct extract_plan));

        ep->data = NULL;

        RUN(galloc(&ep->data,100000,10));

        int len1 , len2;

        RUN(get_dim1(ep->data,&len1));
        RUN(get_dim2(ep->data,&len2));
        int c = 0;
        for(i = 0; i< len1;i++){
                for(j = 0; j < len2;j++){
                        ep->data[i][j] = c;
                        c++;
                }
        }
        /* create a HDF5 file */
        file_id = H5Fcreate ("ex_lite1.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        dims[0] = len1;
        dims[1] = len2;

        /* create and write an integer type dataset named "dset" */
        H5LTmake_dataset(file_id,"/dset",2,dims,H5T_NATIVE_INT,&ep->data[0][0]);

        /* close file */
        H5Fclose (file_id);

        gfree(ep->data);
        MFREE(ep);
        return EXIT_SUCCESS;
ERROR:
        return EXIT_FAILURE;

}
