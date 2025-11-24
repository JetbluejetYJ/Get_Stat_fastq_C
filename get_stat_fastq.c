#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <dirent.h>
#include <sys/stat.h>
#include <zlib.h>
#include "kseq.h"

#define MAX_PATH 1024
#define ZLIB_BUFSIZE (8*1024*1024) // 8MB

typedef struct {
    char sample_prefix[MAX_PATH];
    char fq1[MAX_PATH];
    char fq2[MAX_PATH];
    char out[MAX_PATH];
} ThreadArg;

typedef struct {
    long totalBase, totalRead, base30Count, base20Count, a, c, g, t, n;
    long m1_totalBase, m1_totalRead, m1_base30Count, m1_base20Count, m1_a, m1_c, m1_g, m1_t, m1_n;
    long m2_totalBase, m2_totalRead, m2_base30Count, m2_base20Count, m2_a, m2_c, m2_g, m2_t, m2_n;
} StatResult;

KSEQ_INIT(gzFile, gzread)

void* process_fastq(void* argptr);
void write_output(const char* sample_prefix, const char* out, StatResult* st_result, int paired);

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: %s [path of order data]\n", argv[0]);
        return 1;
    }
    char* dir = argv[1];
    DIR* dp = opendir(dir);
    if (!dp) { perror("opendir"); return 1; }

    struct dirent* ep;
    pthread_t threads[128];
    ThreadArg args[128];
    int tcnt = 0;

    while ((ep = readdir(dp))) {
        if (strstr(ep->d_name, "_1.fastq.gz")) {
            char sample_prefix[MAX_PATH];
            strncpy(sample_prefix, ep->d_name, strlen(ep->d_name) - strlen("_1.fastq.gz"));
            sample_prefix[strlen(ep->d_name) - strlen("_1.fastq.gz")] = 0;

            snprintf(args[tcnt].sample_prefix, MAX_PATH, "%s/%s", dir, sample_prefix);
            snprintf(args[tcnt].fq1, MAX_PATH, "%s/%s_1.fastq.gz", dir, sample_prefix);
            snprintf(args[tcnt].fq2, MAX_PATH, "%s/%s_2.fastq.gz", dir, sample_prefix);
            snprintf(args[tcnt].out, MAX_PATH, "%s/%s.sqs", dir, sample_prefix);

            struct stat st;
            if (stat(args[tcnt].out, &st) == 0 && st.st_size > 0) {
                printf("%s already exists.\n", args[tcnt].out);
                continue;
            }

            pthread_create(&threads[tcnt], NULL, process_fastq, &args[tcnt]);
            tcnt++;
        }
    }
    closedir(dp);

    for (int i = 0; i < tcnt; i++) pthread_join(threads[i], NULL);
    return 0;
}

void* process_fastq(void* argptr) {
    ThreadArg* arg = (ThreadArg*)argptr;
    StatResult st_result = {0};
    int fqScale = 33;
    int paired = 0;

    struct stat st;
    if (stat(arg->fq2, &st) == 0) paired = 1;

    gzFile fp1 = gzopen(arg->fq1, "rb");
    if (!fp1) pthread_exit(NULL);
    gzbuffer(fp1, ZLIB_BUFSIZE); // zlib 버퍼 확장
    kseq_t *ks1 = kseq_init(fp1);

    gzFile fp2 = NULL;
    kseq_t *ks2 = NULL;
    if (paired) {
        fp2 = gzopen(arg->fq2, "rb");
        if (!fp2) pthread_exit(NULL);
        gzbuffer(fp2, ZLIB_BUFSIZE);
        ks2 = kseq_init(fp2);
    }

    while (1) {
        int r1 = kseq_read(ks1);
        int r2 = paired ? kseq_read(ks2) : 0;
        if (r1 < 0 || (paired && r2 < 0)) break;

        // 전체
        st_result.totalBase += ks1->seq.l + (paired ? ks2->seq.l : 0);
        st_result.totalRead += paired ? 2 : 1;
        // R1
        st_result.m1_totalBase += ks1->seq.l;
        st_result.m1_totalRead += 1;
        // R2
        if (paired) { st_result.m2_totalBase += ks2->seq.l; st_result.m2_totalRead += 1; }

        // Base count
        for (int i = 0; i < ks1->seq.l; i++) {
            char base = ks1->seq.s[i];
            if (base == 'A') { st_result.a++; st_result.m1_a++; }
            else if (base == 'C') { st_result.c++; st_result.m1_c++; }
            else if (base == 'G') { st_result.g++; st_result.m1_g++; }
            else if (base == 'T') { st_result.t++; st_result.m1_t++; }
            else if (base == 'N') { st_result.n++; st_result.m1_n++; }
        }
        if (paired) {
            for (int i = 0; i < ks2->seq.l; i++) {
                char base = ks2->seq.s[i];
                if (base == 'A') { st_result.a++; st_result.m2_a++; }
                else if (base == 'C') { st_result.c++; st_result.m2_c++; }
                else if (base == 'G') { st_result.g++; st_result.m2_g++; }
                else if (base == 'T') { st_result.t++; st_result.m2_t++; }
                else if (base == 'N') { st_result.n++; st_result.m2_n++; }
            }
        }

        // Quality
        for (int i = 0; i < ks1->qual.l; i++) {
            int qv = ks1->qual.s[i];
            if (qv >= 30 + fqScale) { st_result.base30Count++; st_result.m1_base30Count++; }
            if (qv >= 20 + fqScale) { st_result.base20Count++; st_result.m1_base20Count++; }
        }
        if (paired) {
            for (int i = 0; i < ks2->qual.l; i++) {
                int qv = ks2->qual.s[i];
                if (qv >= 30 + fqScale) { st_result.base30Count++; st_result.m2_base30Count++; }
                if (qv >= 20 + fqScale) { st_result.base20Count++; st_result.m2_base20Count++; }
            }
        }
    }

    kseq_destroy(ks1); gzclose(fp1);
    if (paired) { kseq_destroy(ks2); gzclose(fp2); }

    write_output(arg->sample_prefix, arg->out, &st_result, paired);
    pthread_exit(NULL);
}

void write_output(const char* sample_prefix, const char* out, StatResult* st_result, int paired) {
    FILE* f = fopen(out, "w");
    if (!f) return;
    #define SAFE_DIV(a,b) ((b) ? ((float)(a)/(b)*100.0f*100.0f/100.0f) : 0.0f)
    #define SAFE_DIV4(a,b) ((b) ? ((float)(a)/(b)*100.0f*10000.0f/10000.0f) : 0.0f)
    #define AVG_READ(a,b) ((b) ? ((float)(a)/(b)*100.0f/100.0f) : 0.0f)

    float f1 = SAFE_DIV(st_result->base30Count, st_result->totalBase);
    float f2 = SAFE_DIV(st_result->base20Count, st_result->totalBase);
    float f3 = SAFE_DIV4(st_result->n, st_result->totalBase);
    float f4 = SAFE_DIV(st_result->g + st_result->c, st_result->totalBase);
    float f5 = AVG_READ(st_result->totalBase, st_result->totalRead);

    fprintf(f, "%s\t%ld\t%ld\t%.4f\t%.2f\t%.2f\t%.2f\n",
        strrchr(sample_prefix, '/') ? strrchr(sample_prefix, '/')+1 : sample_prefix,
        st_result->totalBase, st_result->totalRead, f3, f4, f2, f1);
    fprintf(f, "SampleName : %s\n", strrchr(sample_prefix, '/') ? strrchr(sample_prefix, '/')+1 : sample_prefix);
    fprintf(f, "Total A : %ld\n", st_result->a);
    fprintf(f, "Total C : %ld\n", st_result->c);
    fprintf(f, "Total G : %ld\n", st_result->g);
    fprintf(f, "Total T : %ld\n", st_result->t);
    fprintf(f, "Total N : %ld\n", st_result->n);
    fprintf(f, "Q30 Bases : %ld\n", st_result->base30Count);
    fprintf(f, "Q20 Bases : %ld\n", st_result->base20Count);
    fprintf(f, "Avg.ReadSize : %.2f\n", f5);

    if (paired) {
        f1 = SAFE_DIV(st_result->m1_base30Count, st_result->m1_totalBase);
        f2 = SAFE_DIV(st_result->m1_base20Count, st_result->m1_totalBase);
        f3 = SAFE_DIV4(st_result->m1_n, st_result->m1_totalBase);
        f4 = SAFE_DIV(st_result->m1_g + st_result->m1_c, st_result->m1_totalBase);
        f5 = AVG_READ(st_result->m1_totalBase, st_result->m1_totalRead);
        fprintf(f, "%s_R1\t%ld\t%ld\t%.4f\t%.2f\t%.2f\t%.2f\n",
            strrchr(sample_prefix, '/') ? strrchr(sample_prefix, '/')+1 : sample_prefix,
            st_result->m1_totalBase, st_result->m1_totalRead, f3, f4, f2, f1);
        fprintf(f, "SampleName : %s\n", strrchr(sample_prefix, '/') ? strrchr(sample_prefix, '/')+1 : sample_prefix);
        fprintf(f, "Total A : %ld\n", st_result->m1_a);
        fprintf(f, "Total C : %ld\n", st_result->m1_c);
        fprintf(f, "Total G : %ld\n", st_result->m1_g);
        fprintf(f, "Total T : %ld\n", st_result->m1_t);
        fprintf(f, "Total N : %ld\n", st_result->m1_n);
        fprintf(f, "Q30 Bases : %ld\n", st_result->m1_base30Count);
        fprintf(f, "Q20 Bases : %ld\n", st_result->m1_base20Count);
        fprintf(f, "Avg.ReadSize : %.2f\n", f5);

        f1 = SAFE_DIV(st_result->m2_base30Count, st_result->m2_totalBase);
        f2 = SAFE_DIV(st_result->m2_base20Count, st_result->m2_totalBase);
        f3 = SAFE_DIV4(st_result->m2_n, st_result->m2_totalBase);
        f4 = SAFE_DIV(st_result->m2_g + st_result->m2_c, st_result->m2_totalBase);
        f5 = AVG_READ(st_result->m2_totalBase, st_result->m2_totalRead);
        fprintf(f, "%s_R2\t%ld\t%ld\t%.4f\t%.2f\t%.2f\t%.2f\n",
            strrchr(sample_prefix, '/') ? strrchr(sample_prefix, '/')+1 : sample_prefix,
            st_result->m2_totalBase, st_result->m2_totalRead, f3, f4, f2, f1);
        fprintf(f, "SampleName : %s\n", strrchr(sample_prefix, '/') ? strrchr(sample_prefix, '/')+1 : sample_prefix);
        fprintf(f, "Total A : %ld\n", st_result->m2_a);
        fprintf(f, "Total C : %ld\n", st_result->m2_c);
        fprintf(f, "Total G : %ld\n", st_result->m2_g);
        fprintf(f, "Total T : %ld\n", st_result->m2_t);
        fprintf(f, "Total N : %ld\n", st_result->m2_n);
        fprintf(f, "Q30 Bases : %ld\n", st_result->m2_base30Count);
        fprintf(f, "Q20 Bases : %ld\n", st_result->m2_base20Count);
        fprintf(f, "Avg.ReadSize : %.2f\n", f5);
    }
    fclose(f);
}
