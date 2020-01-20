#ifndef PST_STRUCTS_H
#define PST_STRUCTS_H


/* flat pst structure  */
/* idea : 4*N ints and 4*N floats  */

struct fpst{
        float** prob;
        int** links;
        int l;                  /* length */
        int m;                  /* what is malloced */
};

struct pst {
//        struct pst_node* pst_root;
        //struct pst_node* ppt_root;
        struct fpst* fpst_root;
        struct fpst* fppt_root;
        float p_min;
        float gamma_min;
        float r;
        int L;
};


#endif
