// State-of-the-art kcore algorithm

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>

#include "core/graph.hpp" // gemini system

const VertexId kcore_num = 2;

VertexId compute(Graph<Empty> * graph) {
	MPI_Barrier(MPI_COMM_WORLD);

	double exec_time = - get_time();
	VertexId * degree = graph->alloc_vertex_array<VertexId>();
	VertexSubset * active = graph->alloc_vertex_subset();
	VertexSubset * next_active = graph->alloc_vertex_subset();
	VertexId num_next_actives;
	VertexId num_actives = 0;
	VertexId num_remained_vertices;

	// initialization
	active->clear();
	next_active->clear();
#pragma omp parallel for
	for (VertexId v_i = graph->partition_offset[graph->partition_id]; v_i < graph->partition_offset[graph->partition_id + 1]; ++ v_i) {
		degree[v_i] = graph->out_degree[v_i];
		if (degree[v_i] < kcore_num) {
			active->set_bit(v_i);
			__sync_fetch_and_add(&num_actives, 1);
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &num_actives, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

	while (num_actives) { // if there are vertices removed in the previous iteration, continue the algorithm
		num_next_actives = 0;
		next_active->clear();

		auto sparse_signal = [&](VertexId src) {
			graph->emit(src, 1);
		};
		auto sparse_slot = [&](VertexId src, VertexId msg, VertexAdjList<Empty> outgoing_adj) {
			VertexId actived = 0;
			for (AdjUnit<Empty> * ptr = outgoing_adj.begin; ptr != outgoing_adj.end; ++ ptr) {
				VertexId dst = ptr->neighbour;
				if (__sync_fetch_and_add(&degree[dst], -1) == kcore_num) {
					next_active->set_bit(dst);
					++ actived;
				}
			}
			return actived;

		};
		auto dense_signal = [&](VertexId dst, VertexAdjList<Empty> incoming_adj) {
			VertexId count = 0;
			for (AdjUnit<Empty> * ptr = incoming_adj.begin; ptr != incoming_adj.end; ++ ptr) {
				VertexId src = ptr->neighbour;
				if (active->get_bit(src)) { // ?1: next_active or active // ?2: src or dst
					//
				}
			}
			graph->emit(dst, count); // ?1: src or dst // ?2
		};
		auto dense_slot = [&](VertexId dst, VertexId msg) {
			VertexId current_degree = __sync_fetch_and_add(&degree[dst], -msg);
			if (current_degree >= kcore_num && current_degree - msg < kcore_num) {
				next_active->set_bit(dst); // next_active or active
				return 1;
			}
			return 0;
		};

		num_next_actives = graph->process_edges<VertexId, VertexId>(
				sparse_signal, sparse_slot, dense_signal, dense_slot, active
				);
		if (num_next_actives == 0) {
			break;
		}
		std::swap(active, next_active);
		num_actives = num_next_actives;
	}

	VertexId num_deleted_vertices = 0;
#pragma omp parallel for
	for (VertexId v_i = graph->partition_offset[graph->partition_id]; v_i < graph->partition_offset[graph->partition_id + 1]; ++ v_i) {
		if (degree[v_i] < kcore_num) {
			__sync_fetch_and_add(&num_deleted_vertices, 1);
		}
	}
	VertexId global_num_deleted_vertices;
	MPI_Allreduce(&num_deleted_vertices, &global_num_deleted_vertices, 1, get_mpi_data_type<VertexId>(), MPI_SUM, MPI_COMM_WORLD);
	num_remained_vertices = graph->vertices - global_num_deleted_vertices;

	exec_time += get_time();

	if (! graph->partition_id) {
		printf("exec_time = %.4f(s)\n", exec_time);
		printf("remained vertices = %u\n", num_remained_vertices);
	}

	graph->dealloc_vertex_array(degree);
	delete active;
	delete next_active;
}

int main(int argc, char ** argv) {
	MPI_Instance mpi(&argc, &argv);

	if (argc < 3) {
		printf("kcore_algorithm_1 [file] [vertices]");
		exit(-1);
	}

	Graph<Empty> * graph;
	graph = new Graph<Empty>();
	graph->load_directed(argv[1], std::atoi(argv[2])); // load graph from the file system

	compute(graph);
	int runs = 5;
	for (int i = 0; i < runs; ++ i) { // repeat the experiment 5 times to get the average runtime
		compute(graph);
	}

	delete graph;
	return 0;
}
