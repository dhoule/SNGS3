

#include "dbscan.h"
#include "utils.h"
#include "kdtree2.hpp"

static void usage(char *argv0) {
  const char *params =
    "Usage: %s [switches] -i filename -b -m minpts -n seeds -e epsilon -o output -t threads\n"
    " -i filename : file containing input data to be clustered\n"
    " -b    : input file is in binary format (default no)\n"
    " -m minpts : input parameter of DBSCAN, min points to form a cluster, e.g. 2\n"
    " -n seeds    : input parameter of Sow-n-Grow, number of seeds to be chosen from each partition of data\n"
    " -e epsilon  : input parameter of DBSCAN, radius or threshold on neighbourhoods retrieved, e.g. 0.8\n"
    " -o output : clustering results, format, (each line, point id, clusterid)\n"
    " -t threads  : number of threads to be employed\n\n";
  fprintf(stderr, params, argv0);
  exit(-1);
}

int main(int argc, char** argv) {
  double  seconds;
  int   opt;
  int   minPts, threads;
  int seeds;
  int classical = 0;
  double  eps;
  char*   outfilename = NULL;
  int     isBinaryFile;
  char*   infilename = NULL;
  // some default values
  minPts    = -1;
  eps   = -1;
  seeds = -1;
  isBinaryFile  = 0;
  outfilename   = NULL;
  infilename  = NULL;
  threads   = -1;
  while ((opt=getopt(argc,argv,"i:t:d:p:m:n:e:o:v:z:bxghncul"))!= EOF) {
    switch (opt) {
      case 'i':
        infilename = optarg;
        break;
      case 't':
        threads = atoi(optarg);
        break;
      case 'b':
        isBinaryFile = 1;
        break;
      case 'm':
        minPts = atoi(optarg);
        break;
      case 'n':
        seeds = atoi(optarg);
        break;
      case 'e':
        eps  = atof(optarg);
        break;
      case 'o':
        outfilename = optarg;
        break;
      case 'c':
        classical = 1;
        break;
      case '?':
        usage(argv[0]);
        break;
      default:
        usage(argv[0]);
        break;
    }
  }

  if(infilename == NULL || minPts < 0 || eps < 0 || seeds < 0 || threads < 1) {
    usage(argv[0]);
    exit(-1);
  }

  omp_set_num_threads(threads);
  NWUClustering::ClusteringAlgo dbs;
  double start;
  cout << "\n\nDataset used: " << infilename << endl;
  cout << "Input parameters " << " MinPts " << minPts << " Epsilon " << eps << endl;
  ofstream outputfile;
  outputfile.open("runs.txt");
  
  dbs.set_dbscan_params(eps, minPts, seeds);
  if(dbs.read_file(infilename, isBinaryFile) == -1)
    exit(-1);

  cout << "Dimensions of each point: " << dbs.m_pts->m_i_dims << endl;
  dbs.build_kdtree();
  start = omp_get_wtime();
  run_dbscan_algo_uf(dbs);
  cout << "DBSCAN (total) took " << omp_get_wtime() - start << " seconds.\n" << endl;
  dbs.writeClusters_uf(outputfile);
  outputfile << endl;


  outputfile.close();
  /*if(outfilename != NULL)
  {
    ofstream outfile;
    outfile.open(outfilename);
    dbs.(outfile);
    dbs.writeClusters(outfile);
    outfile.close();
  }*/
  return 0;
}
