/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "adh.h"
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*!
 *  \brief     Initializes and opens an AdH SFILE
 *  \author    Corey Trahan, Ph.D.
 *  \bug       none
 *  \warning   none
 *  \copyright AdH
 *
 * @param[inout] file               (SFILE *)  pointer to an AdH SFILE
 * @param[in]  filebase     (char *) the project root name
 * @param[in]  ext1              (char *) the first filename extension
 * @param[in]  ext2              (char *) the second filename extension
 * @param[in]  suffix          (char *) the filename suffix
 * @param[in]  mode              (char *) read/write mode to open file as
 * @param[in]  stopOnNotFind     (int) a flag to stop AdH if the file is not found
 *
 * \note
 */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void sfile_open(SFILE *file, char *filebase, char *ext1, char *ext2, char *suffix, char *mode, int stopOnNotFind) {
    
	char filename[MAXLINE];        /* the file name */
    sfile_init(file);

    // puts the file name together 
    strcpy(filename, filebase);
    if (ext1 != NULL)   strcat(filename, ext1);
	if (ext2 != NULL)   strcat(filename, ext2);
	if (suffix != NULL) strcat(filename, suffix);
    fprintf(stdout,">> opening filename: %s\n",filename);

    // opens a file if it can and returns the pointer to the file
    file->fp = fopen(filename, mode);
    if(file->fp != NULL) return;

    // otherwise it writes an error message
    if (stopOnNotFind == YES) {
        char s[MAXLINE];
        sprintf(s, "FATAL ERROR :: Can't open file %s mode %s.\n", filename, mode);
        tl_error(s);
	} else {
		fprintf(stdout, "WARNING :: file :: %s is not found.  Continuing on though.\n", filename);
	}
}
