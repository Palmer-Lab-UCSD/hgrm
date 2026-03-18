
#include <io.h>




// textio::TextIO::TextIO(FILE *fileid): fid(fileid) {};
// 
// 
// 
// 
// std::unique_ptr<textio::TextIO> textio::open(const char *filename, 
//         const char *mode) {
// 
//     FILE *fid = fopen(filename, mode);
//     if (ferror(fid)) {
//         fclose(fid);
//         return nullptr;
//     }
// 
//     std::unique_ptr<textio::TextIO> tio = std::make_unique<textio::TextIO>(fid);
// 
//     return std::move(tio);
// }
// 
// 
// textio::STATUS wc(textio::TextIO *tio, textio::FileStats *fs) {
//     if (!fs)
//         return textio::INVALID_ARG_ERROR;
// 
//     size_t nchar = 0;
//     size_t nwords = 0;
//     size_t nlines = 0;
//     size_t nblanklines = 0;
// 
//     size_t word_len = 0;
// 
//     FILE *fid = tio->fid;
// 
//     int c = '\0';
//     while ((c = fgetc(fid)) != EOF) {
// 
//         switch (c) {
//             case '\n':
//                 nlines++;
// 
//                 if (word_len == 0)
//                     nblanklines++;
//                 else {
//                     nwords++;
//                     word_len = 0;
//                 }
//                 break;
//             case ';':
//             case ':':
//             case ',':
//             case '!':
//             case '?':
//             case '(':
//             case ')':
//             case '\"':
//             case '\t':
//             case ' ':
//                 if (word_len == 0)
//                     break;
// 
//                 nwords++;
//                 word_len = 0;
//                 
//                 break;
//             default:
//                 nchar++;
//                 word_len++;
//         }
// 
//     }
// 
//     if (ferror(fid)) {
//         fs = nullptr;
//         return tio->bseek() == 0 ? textio::FERROR : textio::FSEEK_ERROR;
//     }
// 
//     if (feof(fid) == 0) {
//         fs = nullptr;
//         return tio->bseek() == 0 ? textio::FEOF_ERROR : textio::FSEEK_ERROR;
//     }
// 
//     fs->nchar = nchar;
//     fs->nwords = nwords;
//     fs->nlines = nlines;
//     fs->nblanklines = nblanklines;
// 
//     return tio->bseek() ? textio::SUCCESS : textio::FSEEK_ERROR;
// }
// 
// 
// 
// textio::STATUS textio::getline(textio::TextIO *tio, textio::Array<char> *buf) {
//     buf->fill('\0');
// 
//     FILE *fid = buf->fid;
// 
//     int c = 0;
//     while ((c = fgetc(fid)) != EOF) {
// 
//         if (c == '\n') {
//             buf->append('\0');
//             return textio::SUCCESS;
// 
//         buf->append(c)
//     }
// 
//     if (ferror(fid))
//         return textio::FERROR;
// 
//     if (feof(fid) == 0)
//         return textio::FEOF_ERROR;
// 
//     return textio::SUCCESS;
// }
