#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief Reads through a parameter file looking for testcases
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] sm (SMODE_SUPER *)  a pointer to a superModel structure that stores AdH model test case data
 * \note cjt -- there can only be 1 per superModel for now
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void read_bc_TESTCASE(SMODEL_SUPER *sm, FILE *fp) {

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
            sm->testcase = (STESTCASE *) tl_alloc(sizeof(STESTCASE), 1);
            stestcase_init(sm->testcase);

            // ++++++++++++++++++++++++++
            // test cases go here
            // ++++++++++++++++++++++++++
            if (strcmp(token, "SW2") == 0) {
                get_next_token(&token);
                if (strcmp(token, "FLUME") == 0) {
                    sm->testcase->init     =  testcase_sw2_flume_init;
                    sm->testcase->write    =  testcase_sw2_flume_write;
                    sm->testcase->finalize =  testcase_sw2_flume_final;
                }
            }

        }
    }
    rewind(fp);
}

