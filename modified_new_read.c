//Modified version of the EVENTS file reading code written by Dr. Ladislav Subr
//Computes the semi-major axis and eccentricity. Saves it in a file called orb_elements.csv



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#define _G      6.673e-8
#define _pc     3.0856e18
#define _AU     1.49597870e13
#define _M_sun  1.989e33
#define _yr     3.1536e7

#define TTOT    0

#define NIPARAMS  1
#define NFPARAMS  1

#define OPT_HEADER  0x0001  // print only headers
#define OPT_QUIET   0x0002  // do not print headers
#define OPT_VERBOSE 0x0004  // print file position before each record
#define OPT_SYS     0x0008  // print output suitable for SYS files
#define OPT_DR      0x0010  // print relative positions of system components
#define OPT_ZERO    0x0020  // blocks with zero members (headers) are included

struct t_myfile
{
	long lastpos;
	int  offset;
	char *filename;
	FILE *in;
};

#define MAX_TARGET 10

int    g_options;
int    target_id[MAX_TARGET];
int    N_min, N_target;
double mscale, rscale, tscale;
double T_min, T_max;
FILE *outfile;

void V_product(double *V1, double *V2, double *VxV)
{
	VxV[0] = V1[1]*V2[2] - V1[2]*V2[1];    // L_x
	VxV[1] = V1[2]*V2[0] - V1[0]*V2[2];    // L_y
	VxV[2] = V1[0]*V2[1] - V1[1]*V2[0];    // L_z
};

int check_block(struct t_myfile *infile)
{
	int   bs, i;

	if (!(infile->in = fopen(infile->filename, "r"))) return 1;
	
	fread(&bs, 4, 1, infile->in);
	fseek(infile->in, NFPARAMS * 8 + (NIPARAMS + 1) * 4, SEEK_SET);
	fread(&i, 4, 1, infile->in);
	
	infile->offset = 4;
	if (bs != i)
	 {
		fseek(infile->in, NFPARAMS * 8 + (NIPARAMS + 2) * 4, SEEK_SET);
		fread(&i, 4, 1, infile->in);
		infile->offset = 8;
	 };
	
	fseek(infile->in, 0, SEEK_SET);
	
	if (bs != i)
	 {
		fprintf(stderr, "Initial block size mismatch in %s\n", infile->filename);
		return 2;
	 };
	
	return 0;
};

int read_block(struct t_myfile *infile)
{
	int    i, j, k, isys, records, makelog;
	int    *id1, *id2, *nsys;
	double semi, pos;
	double params[NFPARAMS];
	double *r[3], *v[3], *m1, *m2, *H;
	double mu,r_norm,v_norm,E,a;
	double ecc_arr[3],ecc,rv_dot;
	char   inout = '+';

	if (fseek(infile->in, infile->offset, SEEK_CUR) ||
		(fread(&records, 4, 1, infile->in) < 1) ||
		(fread(params, sizeof(double), NFPARAMS, infile->in) < NFPARAMS)) return 1;
	
	if ((!records || (g_options & OPT_VERBOSE)) &&
		(g_options & OPT_HEADER) && (params[TTOT] >= T_min) && (params[TTOT] <= T_max))
		printf("\n# T = %.1f  (0x%08lx)\n", params[TTOT], infile->lastpos);
	if (records)
	 {
		if (records < 0) { records = -records; inout = '-'; };
		if ((g_options & OPT_HEADER) || (params[TTOT] < T_min) || (params[TTOT] > T_max))
			fseek(infile->in, 2 * infile->offset + records * (8 + 9 * sizeof(double)), SEEK_CUR);
		else
		 {
			id1 = (int *)malloc(records * sizeof(int));
			id2 = (int *)malloc(records * sizeof(int));
			m1  = (double *)malloc(records * sizeof(double));
			m2  = (double *)malloc(records * sizeof(double));
			H  = (double *)malloc(records * sizeof(double));
			for (i = 0; i < 3; i++)
			 {
				r[i] = (double *)malloc(records * sizeof(double));
				v[i] = (double *)malloc(records * sizeof(double));
			 };
			nsys = id2;
	
			fseek(infile->in, 2 * infile->offset, SEEK_CUR);
			fread(id1, 4, records, infile->in);
			fread(id2, 4, records, infile->in);
			fread(m1, sizeof(double), records, infile->in);
			fread(m2, sizeof(double), records, infile->in);
			fread(H, sizeof(double), records, infile->in);
			for (i = 0; i < 3; i++) fread(r[i], sizeof(double), records, infile->in);
			for (i = 0; i < 3; i++) fread(v[i], sizeof(double), records, infile->in);
			for (isys = 0; isys < records; isys += nsys[isys])
			 {
				makelog = 0;
//			printf("%2d %2d %4d\n", N_min, nsys[0], target_id);
				if (N_min <= nsys[isys])
				 {
					if (N_target)
					 {
						for (k = 0; k < N_target; k++) if (!makelog)
						 {
							for (j = 1; (j < nsys[isys]) && (id2[isys+j] != target_id[k]); j++);
//				printf("%2d  %2d  %2d  %2d\n", j, nsys[0], id1[j-1], target_id);
							if ((j < nsys[isys]) || (id1[isys+j-1] == target_id[k])) makelog = 1;
						 };
					 }
					else makelog = 1;
				 };
				if (makelog)
				 {
					if (g_options & OPT_SYS)
					 {
						printf("#%d %7d %8.3f %8.3f %8.3f\n", nsys[isys], id1[isys], r[0][isys], r[1][isys], r[2][isys]);
						for (j = isys+1; j < isys+nsys[isys]; j++)
						 {
							semi = -0.5 * (m1[j] + m2[j]) / H[j];
							mu=m1[j]+m2[j];

							r_norm=sqrt(r[0][j]*r[0][j] + r[1][j]*r[1][j] + r[2][j]*r[2][j]);
							v_norm=sqrt(v[0][j]*v[0][j] + v[1][j]*v[1][j] + v[2][j]*v[2][j]);
							E=(v_norm*v_norm/2) - mu/r_norm;
							a=-mu/(2*E);
							
							if (a/semi<1.1 || a/semi>0.9)
							{

								rv_dot=(r[0][j]*v[0][j] + r[1][j]*v[1][j] + r[2][j]*v[2][j]);
								for (i=0;i<3;i++)
								{
									ecc_arr[i]=((((v_norm*v_norm)-mu/r_norm)*r[i][j])-(rv_dot*v[i][j]))/mu;	
								};

								ecc=sqrt(ecc_arr[0]*ecc_arr[0] + ecc_arr[1]*ecc_arr[1] + ecc_arr[2]*ecc_arr[2]);
							}
							else ecc=-1.0;
							fprintf(outfile,"%c,%8.3f,%7d,%7d,%7.12f,%7.12f,%8.10f,%8.10f\n",inout, params[TTOT], id1[j], id2[j], m1[j], m2[j], semi,ecc);


							if (g_options & OPT_DR){
								printf("  %8.8f %8.8f %8.8f\n", r[0][j] , r[1][j], r[2][j]);
								printf("  %8.8f %8.8f %8.8f\n", v[0][j], v[1][j], v[2][j]);}

							else
								printf("\n");
						 };
					 }
					else
					 {
						pos = sqrt(r[0][isys]*r[0][isys] + r[1][isys]*r[1][isys] + r[2][isys]*r[2][isys]);
						printf("%c T = %8.3f %7.3f %6d  %d", inout, params[TTOT], pos, id1[isys], nsys[isys]);
						for (j = 1; j < nsys[isys]; j++) printf(" %5d", id2[isys+j]);
						printf(" %5d\n", id1[isys+j-1]);
					 };
				 };
			 };
	
			free(id1);
			free(id2);
			free(m1);
			free(m2);
			free(H);
			for (i = 0; i < 3; i++) { free(r[i]); free(v[i]); };
		 };
	 }
	else if (g_options & OPT_ZERO)
		fseek(infile->in, 2 * infile->offset, SEEK_CUR);
		
	fseek(infile->in, infile->offset, SEEK_CUR);
	infile->lastpos = ftell(infile->in);
	
	return 0;
};

int main(int argc, char **argv)
{
	int    option;
	struct t_myfile infile;

	N_min = 2;
	g_options = 0;
	N_target = 0;
	mscale = 1.0;
	rscale = 1.0;
	T_min = 0.0;
	T_max = 1e10;
	
	while ((option = getopt(argc, argv, "+N:HqvI:M:R:ST:t:DZ")) != -1) switch (option)
	 {
		case 'D' : g_options |= OPT_DR; break;
		case 'H' : g_options |= OPT_HEADER; break;
		case 'I' : if (N_target < MAX_TARGET) { target_id[N_target] = atol(optarg); N_target++; }; break;
		case 'M' : mscale = atof(optarg); break;
		case 'N' : N_min = atol(optarg); break;
		case 'R' : rscale = atof(optarg); break;
		case 'S' : g_options |= OPT_SYS; break;
		case 'T' : T_max = atof(optarg); break;
		case 'q' : g_options |= OPT_QUIET; break;
		case 't' : T_min = atof(optarg); break;
		case 'v' : g_options |= OPT_VERBOSE; break;
		case 'Z' : g_options |= OPT_ZERO; break;
		case ':' :
		case '?' : return 1;
	 };

	tscale = sqrt(_G * mscale * _M_sun / (rscale * _pc)) / _yr;
	 
	if (argc == optind) { fprintf(stderr, "No files to process. Exitting.\n"); return 1; };

	infile.filename = argv[optind];
	infile.lastpos = 0;
	if (check_block(&infile)) return 1;

	outfile=fopen("orb_elements.csv","w");
	
	printf("\n# %s (%d, %ld)\n", infile.filename, infile.offset, infile.lastpos);
	while(!read_block(&infile));
	fclose(infile.in);
	fclose(outfile);
	return 0;
};
