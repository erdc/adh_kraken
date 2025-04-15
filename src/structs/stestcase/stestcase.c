/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*! \file  stestcase.c Collects routine for AdH models test cases */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Initializes a stestcase structure
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] tc (STESTCASE *)  a pointer to a STESTCASE structure that stores AdH model test case data
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void stestcase_init(STESTCASE *tc) {

    tc->testcase_init  = NULL;
    tc->testcase_write = NULL;
    tc->testcase_final = NULL;
    
} 

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Reads through a parameter file looking for testcases
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] tc (STESTCASE *)  a pointer to a STESTCASE structure that stores AdH model test case data
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void stestcase_read(STESTCASE *tc, FILE *fp) {
    
    int i;
    size_t len = 0;
    ssize_t read;
    char *line = NULL, *token = NULL, str[MAXLINE] = "";
    
    // +++++++++++++++++++++++++++++++++++
    // assign pointers to testcase functions
    // +++++++++++++++++++++++++++++++++++
    rewind(fp);
    while ((read = getline(&line, &len, fp)) != -1) {
        get_token(line,&token); if (token == NULL) continue;
        if (strcmp(token, "TEST") == 0) {
            get_next_token(&token);
            
            // ++++++++++++++++++++++++++
            // test cases go here
            // ++++++++++++++++++++++++++
            if (strcmp(token, "SW2 FLUME") == 0) {
                tc->testcase_init  =  testcase_sw2_flume_init;
                tc->testcase_write =  testcase_sw2_flume_write;
                tc->testcase_final =  testcase_sw2_flume_final;
            }
            
        }
    }
    
    
    rewind(fp);
}
