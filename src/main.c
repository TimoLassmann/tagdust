/*

  Copyright (C) 2013 Timo Lassmann <timolassmann@gmail.com>

  This file is part of TagDust.

  TagDust is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  TagDust is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Tagdust.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "scTagDust.h"
#include "nuc_code.h"


#include "interface.h"
/*#include "nuc_code.h"
#include "io.h"
#include "misc.h"
#include "tagdust2.h"
#include "barcode_hmm.h"
#include <math.h>
*/

int main (int argc,char * argv[]) {
        struct parameters* param = NULL;
        int i;
        RUN(init_nuc_code());

        RUN(interface(&param,argc,argv));

        if(!param){
                return EXIT_SUCCESS;
        }
        if(param->read_structure->num_segments == 0 && param->arch_file == NULL){
                ERROR_MSG("ERROR: No read architecture found.\n");
        }

        RUN(QC_read_structure(param));

        ASSERT(param->infiles > 0, "Number of inputs has to be greater than 0");

        ASSERT(param->outfile != NULL, "No output file suffic");

        if(param->arch_file != NULL){
                ASSERT(my_file_exists(param->arch_file),"Could not find file %s.",param->arch_file);
        }
        if(param->infiles){
                for(i = 0; i < param->infiles;i++){
                        ASSERT(my_file_exists(param->infile[i]),"Could not find file %s.",param->infile[i]);
                }
        }

        //sprintf(param->buffer,"Start Run\n--------------------------------------------------\n");
        //param->messages = append_message(param->messages, param->buffer);
        //hmm_controller_multiple(param);
        free_param(param);

        return EXIT_SUCCESS;
ERROR:
        if(param){
                //fprintf(stdout,"%s",param->errmsg);
                free_param(param);
        }
        return EXIT_SUCCESS;
}












