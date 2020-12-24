
#ifndef _DBSCAN_

#define _DBSCAN_



#include "utils.h"

#include "clusters.h"

namespace NWUClustering
{

	class ClusteringAlgo : public Clusters

	{

	public:

		ClusteringAlgo(){ }

		virtual ~ClusteringAlgo();

		// functions for dbscan algorithm
		void set_dbscan_params(double eps, int minPts, int seeds);

		void 	writeClusters(ostream& o); // regular dbscan algorithm

		void    writeClusters_uf(ostream& o); // union find dbscan algorithm

	public:

		// parameters to run dbscan algorithm

		double 	m_epsSquare;

		int 	m_minPts;

		int m_seeds;

		//int     m_parcent_of_data;
		// noise vector

        	vector<bool> m_noise;

	       	// noise vector

        	vector<bool> m_visited;



		vector <int> m_parents;

		vector <int> m_corepoint;

		vector <int> m_member;

	};	



	void run_dbscan_algo_uf(ClusteringAlgo& dbs); // union find dbscan algorithm

	void run_dbscan_algo(ClusteringAlgo& dbs); // regular dbscan algorithm

};



#endif
