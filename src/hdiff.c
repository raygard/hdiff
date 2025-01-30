// hdiff.c -- histogram diff
// Copyright 2025 Ray Gardner
// License: 0BSD
// vi: tabstop=2 softtabstop=2 shiftwidth=2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>

#include <unistd.h>
#include <sys/stat.h>

#ifndef minof
#define minof(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef maxof
#define maxof(a, b) ((a) > (b) ? (a) : (b))
#endif

static void error_exit(char *s)
{
  fprintf(stderr, "%s\n", s);
  exit(2);
}

static void xstat(char *path, struct stat *buf)
{
  if (stat(path, buf)) fprintf(stderr, "can't stat %s\n", path);
}

static void *xzalloc(size_t n)
{
  void *p = calloc(1, n);
  if (!p) error_exit("out of memory");
  return p;
}

static void *xrealloc(void *p, size_t n)
{
  p = realloc(p, n);
  if (!p) error_exit("out of memory");
  return p;
}

static char *read_input(FILE *fp, size_t *plen)
{
  size_t bufsize = 4096, k = 0;
  long len = 0;
  char *buf;
  if (!fseek(fp, 0, SEEK_END)) {
    // Seekable file.
    if ((len = ftell(fp)) < 0 || len >= INT_MAX) return 0;
    rewind(fp);
    buf = xzalloc(len + 1);
    if ((long)fread(buf, 1, len, fp) != len) {
      free(buf);
      return 0;
    }
    buf[*plen = len] = 0;
    return buf;
  }
  // Non seekable; e.g. stdin from terminal or pipe
  // Need to read into expandable buffer until EOF
  buf = xzalloc(bufsize + 1);
  for (;;) {
    k = fread(buf + len, 1, bufsize - len, fp);
    len += k;
    if (len < (long)bufsize) break;
    // buffer full
    bufsize = bufsize / 5 * 8;
    buf = xrealloc(buf, bufsize + 1);
  }
  if (ferror(fp)) error_exit("i/o error");
  buf[*plen = len] = 0;
  return buf;
}

static int eqln(char *a, int alen, char *b, int blen)
{
  return alen == blen && !memcmp(a, b, alen);
}

static int fill_offset_vec(int *v, char *dat, size_t len)
{
  int k = 1;
  char *p, *b = dat;
  while ((p = memchr(b, '\n', len))) {  // SAFE?
    len -= ++p - b;
    b = p;
    k++;    // one more than num of newlines seen
    if (v) v[k] = p - dat;  // v[2] == offset to line 2 etc.
  }
  if (len) v ? (v[++k] = b - dat + len) : k++;
  return k - 1;
}

static unsigned hash(char *s, int len)
{
  unsigned h = 5381;    // djb2a hash
  // while (len--) h = h * 33 ^ (unsigned char)*s++;
  while (len--) // ignore space and case
    if (isspace((unsigned char)*s)) s++;
    else h = h * 33 ^ tolower((unsigned char)*s++);
  return h;
}

static int hash_find(char *s, int len,
        int *hashtbl, unsigned hmask, char *dat, int *offs)
{
  int k;
  unsigned h, h0 = hash(s, len);
  h = h0 & hmask;
  while ((k = hashtbl[h]) && !eqln(s, len, dat + offs[k], offs[k+1] - offs[k]))
    h = (h * 5 + 1 + (h0 >>=5)) & hmask;
  // Now hashtbl[h] is 0 or it's a hit
  return h;
}

static void push_quad(int **stk, int *stkmax, int *stkcnt,
        int a, int b, int c, int d)
{
  if (*stkcnt + 4 > *stkmax)
    *stk = xrealloc(*stk, sizeof(**stk)*(*stkmax = *stkmax * 8 / 5));
  (*stk)[(*stkcnt)++] = a;
  (*stk)[(*stkcnt)++] = b;
  (*stk)[(*stkcnt)++] = c;
  (*stk)[(*stkcnt)++] = d;
}

static void pop_quad(int *stk, int *stkcnt, int *a, int *b, int *c, int *d)
{
  if (*stkcnt < 4) error_exit("stack underflow");
  *d = stk[--*stkcnt];
  *c = stk[--*stkcnt];
  *b = stk[--*stkcnt];
  *a = stk[--*stkcnt];
}

// Formats time_t in ISO format ("yyyy-mm-dd hh:mm:ss.nnnnnnnnn -hhmm").
static char *fmt_iso_time(char *buf, size_t len, time_t *mtime, long nsec)
{
  char *s = buf;
  s += strftime(s, len, "%Y-%m-%d %H:%M:%S", localtime(mtime));
  sprintf(s, ".%9.9ld", nsec);
  s += 10;
  s += strftime(s, len, " %z", localtime(mtime));
  return buf;
}

static void printhdr(char *prefix, char *fn)
{
  char tm[64];
  struct stat statbuf;
  struct timespec ts;
  time_t mtim;
  long nsec;
  if (strcmp(fn, "-")) {
    xstat(fn, &statbuf);
    mtim = statbuf.st_mtime;
    nsec = statbuf.st_mtim.tv_nsec;
  } else {
    mtim = time(NULL);
    clock_gettime(CLOCK_REALTIME, &ts);
    nsec = ts.tv_nsec;
  }
  printf("%s %s\t%s\n", prefix, fn, fmt_iso_time(tm, sizeof(tm), &mtim, nsec));
}

static void slider(int *d, int j, int k, int *offs, char *dat, int nrecs)
{
  int lenlo, lenhi, n;
  while (d[k+1] < nrecs && d[k+1] < d[k+4]) {
    lenlo = offs[d[k]+1]-offs[d[k]];
    lenhi = offs[d[k+1]+1]-offs[d[k+1]];
    if (!eqln(dat+offs[d[k]], lenlo, dat+offs[d[k+1]], lenhi)) break;
    for (n = j; n < j + 4; n++) d[n]++;
  }
}

static void print_line(int flg, char *dat, int *offs, int i)
{
  int k = offs[i+1] - offs[i];
  printf("%c %*.*s", flg, k, k, dat+offs[i]);
}

static void print_diff(char *kind, char *adat, int *aoffs,
        char *bdat, int *boffs, int alo, int ahi, int blo, int bhi)
{
  char dbuf[24], abuf[24];
  if (alo>=--ahi) snprintf(dbuf, sizeof dbuf, "%d", ahi);
  else snprintf(dbuf, sizeof dbuf, "%d,%d", alo, ahi);
  if (blo>=--bhi) snprintf(abuf, sizeof abuf, "%d", bhi);
  else snprintf(abuf, sizeof abuf, "%d,%d", blo, bhi);
  printf("%s%s%s\n", dbuf, kind, abuf);
  while (alo <= ahi)
    print_line('<', adat, aoffs, alo++);
  if (*kind == 'c') printf("---\n");
  while (blo <= bhi)
    print_line('>', bdat, boffs, blo++);
}

static void print_diffs_default(int *d, int diffscnt,
        char *adat, int *aoffs, int arecs, char *bdat, int *boffs, int brecs)
{
  int j;
  for (j = 4; j < diffscnt-4; j += 4) {
    if (d[j] == d[j+1]) slider(d, j, j+2, boffs, bdat, brecs);
    if (d[j+2] == d[j+3]) slider(d, j, j, aoffs, adat, arecs);
  }
  for (j = 4; j < diffscnt-4; j += 4) {
    if (d[j] < d[j+1])
      if (d[j+2] < d[j+3])
        print_diff("c", adat, aoffs, bdat, boffs, d[j], d[j+1], d[j+2], d[j+3]);
      else
        print_diff("d", adat, aoffs, bdat, boffs, d[j], d[j+1], d[j+2], d[j+3]);
    else if (d[j+2] < d[j+3])
      print_diff("a", adat, aoffs, bdat, boffs, d[j], d[j+1], d[j+2], d[j+3]);
    else
      printf("NO DATA IN CHANGE? %d %d %d %d\n", d[j], d[j+1], d[j+2], d[j+3]);
  }
}

static void print_line_u(char *prefix, int i, char *dat, int *offs)
{
  int k = offs[i+1] - offs[i];
  printf("%s%*.*s", prefix, k, k, dat+offs[i]);
}
static void print_diffs_unified(int *d, int diffscnt,
      char *afn, char *adat, int *aoffs, int arecs,
        char *bfn, char *bdat, int *boffs, int brecs, int nctx)
{
  int i = 4, j, k, n, calo, cahi, cblo, cbhi;
  if (diffscnt < 12) return;
  printhdr("---", afn);
  printhdr("+++", bfn);

  for (j = 4; j < diffscnt-4; j += 4) {
    if (d[j] == d[j+1]) slider(d, j, j+2, boffs, bdat, brecs);
    if (d[j+2] == d[j+3]) slider(d, j, j, aoffs, adat, arecs);
  }
  while (i < diffscnt-4) {
    // Find chunks separated by <= 2 * nctx
    for (j = i + 4; j < diffscnt - 4 && d[j] - d[j-3] <= 2 * nctx; j += 4)
      ;
    calo = maxof(d[i] - nctx, 1);
    cahi = minof(d[j-3] + nctx, d[j]);
    cblo = maxof(d[i+2] - nctx, 1);
    cbhi = minof(d[j-1] + nctx, d[j+3]);

    // Adjust begin line numbers if count is zero
    printf("@@ -%d", cahi-calo?calo:calo-1);
    if(cahi-calo!=1)printf(",%d", cahi-calo);

    printf(" +%d", cbhi-cblo?cblo:cblo-1);
    if(cbhi-cblo!=1)printf(",%d", cbhi-cblo);
    printf(" @@\n");

    for (k = i; k < j; k+=4) {
      for (n=(k>i?d[k-3]:maxof(1,d[k]-nctx)); n < d[k]; n++)
        print_line_u(" ", n, adat, aoffs);
      for (n = d[k]; n < d[k+1]; n++) print_line_u("-", n, adat, aoffs);
      for (n = d[k+2]; n < d[k+3]; n++) print_line_u("+", n, bdat, boffs);
    }
    for (n = d[j-3]; n < minof(d[j-3]+nctx, d[j]); n++)
      print_line_u(" ", n, adat, aoffs);
    i = j;
  }
}

static void print_diffs(int *d, int diffscnt,
      char *afn, char *adat, int *aoffs, int arecs,
        char *bfn, char *bdat, int *boffs, int brecs, int dtype, int nctx)
{
    if (dtype == 'u') print_diffs_unified(d, diffscnt,
            afn, adat, aoffs, arecs, bfn, bdat, boffs, brecs, nctx);
    else print_diffs_default(d, diffscnt,
            adat, aoffs, arecs, bdat, boffs, brecs);
}

static int getcnt(int ap, int *acnt)
{
  return acnt[ap] >= 0 ? acnt[ap] : acnt[-acnt[ap]];
}

static int beqa(int ap, int *acnt, int bp, int *bref)
{
  return bref[bp] == (acnt[ap] >= 0 ? ap : -acnt[ap]);
}

#define INITLOWCNT 512      // jgit uses 65
static int find_best_matching_region(int *acnt, int *anext, int *bref,
        int alo, int ahi, int blo, int bhi,
        int *kalo, int *kahi, int *kblo, int *kbhi)
{
  int i, j, lowcnt, rgnlowcnt, cnt, nextj, nexti, nalo, nahi, nblo, nbhi;
  *kalo = 0;
  if (alo == ahi || blo == bhi) return 0;
  lowcnt = INITLOWCNT;
  for (i = blo; i < bhi; i = nexti) {
    nexti = i + 1;
    if (!(j = bref[i])) continue;
    if (j >= ahi || !acnt[j] || acnt[j] > lowcnt) continue;
    // try to find j inside A range
    while (j < alo && anext[j]) j = anext[j];
    if (j < alo || j >= ahi) continue;
    for (;;) {
      nextj = anext[j];
      // alo <= j < ahi and a[j] matches b[i]
      // set current match, then expand it
      nalo = j; nahi = j+1; nblo = i; nbhi = i+1;
      rgnlowcnt = getcnt(nalo, acnt);
      while (alo < nalo && blo < nblo && beqa(nalo-1, acnt, nblo-1, bref)) {
        nalo--; nblo--;
        cnt = getcnt(nalo, acnt);
        if (cnt < rgnlowcnt) rgnlowcnt = cnt;
      }
      while (nahi < ahi && nbhi < bhi && beqa(nahi, acnt, nbhi, bref)) {
        cnt = getcnt(nahi, acnt);
        if (cnt < rgnlowcnt) rgnlowcnt = cnt;
        nahi++; nbhi++;
      }
      if (!*kalo || nahi - nalo > *kahi - *kalo || rgnlowcnt < lowcnt) {
        *kalo = nalo; *kahi = nahi; *kblo = nblo; *kbhi = nbhi;
        lowcnt = rgnlowcnt;
      }
      if (nexti < nbhi) nexti = nbhi;
      while (nextj && nextj < nahi) nextj = anext[nextj];
      if (!nextj || nextj >= ahi) break;
      j = nextj;
    }
  }
  if (*kalo) return 1;
  return 0;
}

static int diff(char *afn, char *bfn, int dtype, int nctx)
{
  FILE *afp = stdin, *bfp = stdin;
  char *adat, *bdat;    // data areas for files A and B
  size_t alen, blen;    // length of files A and B
  int arecs, brecs;     // rec (line) count for A and B
  int *aoffs, *boffs;   // offset vectors to recs (lines) in data areas
  int alo, ahi, blo, bhi; // bounds of current region
  int kalo, kahi, kblo, kbhi; // bounds of best matching region
  // See Implementation Notes at end
  int *anext, *bref, *acnt;
  int *hashtbl, hashtblsize, hx;
  int rstackcnt = 0, rstackmax = 100; // rstackmax s/b multiple of 4
  int *rstack = xzalloc(rstackmax * sizeof(*rstack));
  int diffscnt = 0, diffsmax = 100; // diffsmax s/b multiple of 4
  int *diffs = xzalloc(diffsmax * sizeof(*diffs));
  int i, j, k;

  // open files A and B
  if (strcmp(afn, "-"))
    afp = fopen(afn, "rb");
  if (strcmp(bfn, "-"))
    bfp = fopen(bfn, "rb");
  if (!afp) error_exit("Can't open file1");
  if (!bfp) error_exit("Can't open file2");
  if (afp == bfp) error_exit("Both inputs cannot be stdin");

  // read files to areas adat and bdat
  adat = read_input(afp, &alen);
  if (!adat) error_exit("cannot read file1");
  if (afp != stdin) fclose(afp);

  bdat = read_input(bfp, &blen);
  if (!bdat) error_exit("cannot read file2");
  if (bfp != stdin) fclose(bfp);

  // first call to fill_offset_vec just counts lines
  arecs = fill_offset_vec(0, adat, alen);
  brecs = fill_offset_vec(0, bdat, blen);

  // alloc vectors and hash table
  aoffs = xzalloc((arecs + 2) * sizeof(*aoffs));
  boffs = xzalloc((brecs + 2) * sizeof(*boffs));
  anext = xzalloc((arecs + 1) * sizeof(*anext));
  bref = xzalloc((brecs + 1) * sizeof(*bref));
  acnt = xzalloc((arecs + 1) * sizeof(*acnt));
  for (hashtblsize = 128; hashtblsize < arecs / 4 * 5; hashtblsize *= 2) {}
  hashtbl = xzalloc(hashtblsize * sizeof(*hashtbl));

  // fill offset vectors pointing to lines in adat and bdat
  arecs = fill_offset_vec(aoffs, adat, alen);
  brecs = fill_offset_vec(boffs, bdat, blen);

  // fill anext and bref vectors using hash table
  // see Implementation Notes at end
  for (k = arecs; k; k--) {   // fill anext[]
    hx = hash_find(adat + aoffs[k], aoffs[k+1] - aoffs[k],
                  hashtbl, hashtblsize - 1, adat, aoffs);
    if (hashtbl[hx]) anext[k] = hashtbl[hx];
    hashtbl[hx] = k;
  }
  for (k = 1; k <= brecs; k++) {    // fill bref[]
    hx = hash_find(bdat + boffs[k], boffs[k+1] - boffs[k],
                  hashtbl, hashtblsize - 1, adat, aoffs);
    bref[k] = hashtbl[hx];
  }
  free(hashtbl);  // we're done with hash table

  // Fill count of A line occurrences vec (acnt)
  for (i = 1; i <= arecs; i++)
    if (!acnt[i])
      for (j = i, acnt[i] = 1; (j = anext[j]); acnt[i]++) acnt[j] = -i;

  push_quad(&diffs, &diffsmax, &diffscnt, 1, 1, 1, 1);
  if (arecs || brecs) // don't push if both files empty!
    push_quad(&rstack, &rstackmax, &rstackcnt, 1, arecs+1, 1, brecs+1);
  while (rstackcnt) {
    pop_quad(rstack, &rstackcnt, &alo, &ahi, &blo, &bhi);
    // adjust acnt[x] for x in A range (alo <= x < ahi)
    for (i = blo; i < bhi; i++)
      if (bref[i]) acnt[bref[i]] = 0;
    for (i = alo; i < ahi; i++)
      if (acnt[i] < 0) acnt[-acnt[i]] = 0;
      else acnt[i] = 0;
    for (i = alo; i < ahi; i++)
      if (acnt[i] < 0) acnt[-acnt[i]]++;
      else if (!acnt[i]) acnt[i]++;
      else error_exit("refill error");
    if (find_best_matching_region(acnt, anext, bref, alo, ahi, blo, bhi,
                &kalo, &kahi, &kblo, &kbhi)) {
      // k-region is a match (kalo upto kahi matches kblo upto kbhi)
      // if region after includes at least 1 line, push
      if (kahi < ahi || kbhi < bhi)
        push_quad(&rstack, &rstackmax, &rstackcnt, kahi, ahi, kbhi, bhi);
      // if region before includes at least 1 line, push
      if (alo < kalo || blo < kblo)
        push_quad(&rstack, &rstackmax, &rstackcnt, alo, kalo, blo, kblo);
    } else {
      // Only non-empty regions are pushed, so this will be nonempty
      push_quad(&diffs, &diffsmax, &diffscnt, alo, ahi, blo, bhi);
    }
  }
  // facilitate printing diffs with this empty final diff region
  push_quad(&diffs, &diffsmax, &diffscnt, arecs+1, arecs+1, brecs+1, brecs+1);
  print_diffs(diffs, diffscnt, afn, adat, aoffs, arecs,
          bfn, bdat, boffs, brecs, dtype, nctx);
  free(adat);
  free(bdat);
  free(aoffs);
  free(boffs);
  free(anext);
  free(bref);
  free(acnt);
  free(diffs);
  free(rstack);
  return diffscnt > 8;
}

static char *version = "25.01 20250129";

int main(int argc, char **argv)
{
  int opt, dtype = 'd', nctx = 3;
  char *usage = "Usage: diff [-u | -U n] file1 file2\n-V version info";

  while ((opt = getopt(argc, argv, "uU:V")) != -1) {
    switch (opt) {
      case 'u':
        dtype = 'u';
        break;
      case 'U':
        dtype = 'u';
        nctx = atoi(optarg);
        break;
      case 'V':
        printf("version (%s), compiled %s %s\n",
                    version, __DATE__, __TIME__);
        exit(2);
        break;
      default:
        error_exit(usage);
    }
  }
  if (optind + 2 != argc) error_exit(usage);

  return diff(argv[optind], argv[optind + 1], dtype, nctx);
}

/* Implementation Notes
 *
 * This is an implementation of "histogram diff" based on the original
 * implementation in jgit. It may differ slightly in output from the
 * original or from 'git diff --histogram'.
 *
 * The basic idea is to find a region of consecutive matching lines
 * between files A and B, then repeat the process recursively on the
 * regions above and below the matching region until encountering
 * regions that contain no matching lines; these are the differences.
 *
 * A region is a set of consecutive lines from file A and file B, and a
 * matching region is one where there are the same number of lines from
 * A and B and they all match in order. A region is held in four
 * variables, e.g. alo, ahi, blo, bhi, where alo is the first line
 * number from A and ahi is the line number *following* the last line of
 * the A side of the region (a half-open range), and similarly for blo
 * and bhi.
 * 
 * Begin by reading the files into data areas adat and bdat.
 * Create vectors (arrays) aoffs and boffs of offsets of line starting
 * positions in the data areas. The fill_offset_vec() function runs
 * twice for each area; the first to count lines (arecs, brecs) and then
 * to fill the offset vectors.
 * Note that lines are numbered from 1 to arec (brec), and position 0
 * of aoffs (boffs) is unused. The aoffs[arecs+1] (boffs[brecs+1])
 * slot will have the file length alen (blen), which is also the
 * offset one beyond the data in adat (bdat).
 *
 * A vector anext[] contains, for each line in A, the line number of the
 * next line in A that matches the current line, or 0 if there is no
 * next matching line. The A lines are traversed from last to first so
 * that the links can go forward. A hash table is used to keep the line
 * number of the most recently encountered occurrence of each line.
 * 
 * A vector bref[] contains, for each line in B, the line number of the
 * first occurrence of the matching line in A, or 0 if no A line
 * matches. The hash table is used to fill bref[], then free()'d.
 * 
 * A vector acnt[] contains, for each line in A, the number of times
 * that line occurs in A (at the first location that line occurs) or a
 * negative index of that first location. That is, if line i is the
 * first occurrence of a line, acnt[i] holds the number of occurrences
 * of that line in A, and if it's a subsequence occurrence of the same
 * line, then acnt[-acnt[i]] holds the number of occurrences.
 * 
 * Two dynamically grown lists are maintained, rstack[] and diffs[].
 * rstack[] is the recursion stack to hold regions that have not yet
 * been examined to try to find a matching region inside. diffs[] will
 * hold the bounds of the difference regions as they are located.
 * 
 * To begin the search for differences, a "dummy" difference is pushed
 * on the diffs[] list to facilitate the writing of the actual diff
 * output, and all of each file (1, arecs+1, 1, brecs+1) is pushed on
 * rstack[].
 * 
 * Then while rstack[] is nonempty, pop a region from it into (alo, ahi,
 * blo, bhi) and try to find a "best" matching region inside that
 * region. If such a match is found, push the region after it, and then
 * the one before it, onto the stack. If no matching region is found,
 * then push (append) the current region onto the diffs[] list. Then
 * continue the loop until the stack is empty.
 *
 * Before searching for the best matching region, adjust the acnt[]
 * values for the current region (alo, ahi...) to reflect just the
 * occurrences of each A line within that region.
 * 
 * After the stack is empty, process the diffs[] list into the actual
 * diff output.
 * 
 * Function find_best_matching_region() is the heart of the algorithm.
 * It takes (alo, ahi, blo, bhi) as the current region as input and
 * returns (kalo, kahi, kblo, kbhi) as the best match, if any. Set
 * lowcount to an initial value high enough to avoid matching lines in A
 * that occur "too many" times (to avoid worse-tha-quadratic time
 * performance on pathological inputs).
 *
 * Traverse the B lines of the current region using index i. If there is
 * a matching line in A, set j to that line (j = bref[i]). If it falls
 * beyond the end of the A region, or if it occurs more than lowcount
 * times in the A region, go to the next line in B and repeat.
 * 
 * If j falls before alo, then use anext[] to find the first A line in
 * the region that matches the B line (i). If none is found, continue
 * with the next line in B.
 * 
 * For each A line (j) in the current region that matches the current B
 * line (i), try to find a "best" match. First, create a region of just
 * the one line (nalo, nahi, nblo, nbhi) set from (j, j+1, i, i+1). Then
 * widen the region looking for matching lines before and after it. Use
 * the acnt[] and bref[] vectors to check for line equality. Also keep
 * track of the lowest count (from acnt[]) for any A line in the
 * matching region.
 * 
 * If this is the first matching region found, or if it is longer than
 * the best matching region found so far, or if the count is lower than
 * the count of the best matching region so far (lowcount), set this
 * region as the new best: set (kalo, kahi, kblo, kbhi) from (nalo,
 * nahi, nblo, nbhi) and set lowcount to the lowest A count in the
 * matching region.
 * 
 * Also, when a match is found, set the next line in B to be examined to
 * nbhi (skipping lines before and within the matching region). But also
 * continue looking at matching lines in A (via anext[]) if any, before
 * going on to another line in B.
 * 
 * After reaching the end of B lines (i == bhi), return 1 if any
 * matching region was found, or 0 to signify no matching region.
 * 
 * Before writing any diff output, jgit "adjusts" some diffs that are
 * insertions or deletions, attempting to make them look better. jgit
 * calls this "normalize"; my implementation calls it slider(). You may
 * have seen diffs where an insertion looks like:
 * 
 * + }
 * +
 * + static int foo()
 * + {
 * +   whatever...
 *   }
 *  
 * + more stuff
 * 
 * The actual change was
 * 
 *   }
 *  
 * + static int foo()
 * + {
 * +   whatever...
 * + }
 * +
 * + more stuff
 * 
 * The slider() function makes this adjustment. It usually makes the
 * diff more intuitive to read, but may not always.
 */
