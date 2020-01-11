#ifndef BASE_QUALITY_H
#define BASE_QUALITY_H

#include <stdint.h>

struct qsubscore;

extern int calc_score_matrix(struct qsubscore** mat,  float base_error, float indel_freq);
extern float get_qsubscore(struct qsubscore* subm, uint8_t a, uint8_t b, uint8_t q);


#endif
