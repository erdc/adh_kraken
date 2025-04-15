#ifndef H_STESTCASE_
#define H_STESTCASE_

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

typedef struct {

    void (*testcase_init)(SMODEL_SUPER *);
    void (*testcase_write)(SMODEL_SUPER *);
    void (*testcase_final)(SMODEL_SUPER *);

} STESTCASE;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// methods

void stestcase_init(STESTCASE *tc);
void stestcase_read(STESTCASE *tc, FILE *fp);


// flume
void testcase_sw2_flume_init(SMODEL_SUPER *sm);
void testcase_sw2_flume_final(SMODEL_SUPER *sm);
void testcase_sw2_flume_write(SMODEL_SUPER *sm);

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif

