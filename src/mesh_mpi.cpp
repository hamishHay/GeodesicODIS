#include "mesh.h"
#include "globals.h"
#include "math.h"
#include "array1d.h"
#include "array2d.h"
#include "array3d.h"
#include "mathRoutines.h"
#include "gridConstants.h"
// #include "sphericalHarmonics.h"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include <filesystem>
#include <vector>


#ifdef _MPI
    #include "mpi.h"

    namespace 
    {
        
        constexpr int TAG_NODE_COUNT   = 10;
        constexpr int TAG_NODE_IDS     = 11;
        constexpr int TAG_FACE_COUNT   = 20;
        constexpr int TAG_FACE_IDS     = 21;
        constexpr int TAG_VERTEX_COUNT = 30;
        constexpr int TAG_VERTEX_IDS   = 31;

        constexpr int TAG_BASE_NODE_DATA    = 100;
        constexpr int TAG_BASE_FACE_DATA    = 200;
        constexpr int TAG_BASE_VERTEX_DATA  = 300;
        constexpr int TAG_BASE_NODE_XYZ_DATA= 400;
    

        struct MPIElementConnectivity {
            const char * name;              // for debug output, e.g. "NODE", "FACE", "VERTEX"

            unsigned num_ng;                // index where ghosts start
            unsigned num;                   // total count (owned + ghost)

            Array2D<unsigned> & region_ID;  // (region, remote_local_id) per entity

            std::vector<unsigned> & recv_regions;
            std::vector<unsigned> & send_regions;

            std::vector<std::vector<unsigned>> & receive_IDs;
            std::vector<std::vector<unsigned>> & send_IDs;
            std::vector<std::vector<unsigned>> & send_IDs_ordered;
            std::vector<std::vector<unsigned>> & send_IDs_map;

            int tag_count;
            int tag_ids;
        };

        // Builds send/receive ID lists and region counts for an element type
        // (nodes, faces, or vertices) for each process.
        void FindMPIElementConnectivity(MPIElementConnectivity & element, int proc_id, unsigned nproc,
                                    std::ostringstream & outstring)
        {
            // std::vector<unsigned> ghost_num_in_region(nproc, 0);

            for (unsigned k = 0; k < nproc; k++) element.recv_regions[k] = 0;

            for (unsigned i = element.num_ng; i < element.num; i++) {
                unsigned k = element.region_ID(i, 0);
                // ghost_num_in_region[k] += 1;
                element.recv_regions[k] += 1;
            }

            for (unsigned k = 0; k < nproc; k++) {
                if (proc_id == (int)k) continue;

                unsigned num_to_receive = element.recv_regions[k];//ghost_num_in_region[k];
                unsigned num_to_send = 0;

                MPI_Status status;
                MPI_Sendrecv(
                    &num_to_receive, 1, MPI_UNSIGNED, k, element.tag_count,
                    &num_to_send,    1, MPI_UNSIGNED, k, element.tag_count,
                    MPI_COMM_WORLD, &status
                );

                element.send_regions[k] = num_to_send;

                outstring << "IN PROCESS " << proc_id << " THERE ARE " << num_to_receive
                        << " GHOST " << element.name << "S FROM PROCESS " << k << std::endl;
                outstring << "PROCESS " << proc_id << " WILL SEND " << num_to_send
                        << " " << element.name << "S TO PROCESS " << k << std::endl;

                std::vector<unsigned> IDs(num_to_receive);
                std::vector<unsigned> IDs_to_request(num_to_receive);
                std::vector<unsigned> IDs_to_send(num_to_send);

                unsigned count = 0;
                for (unsigned i = element.num_ng; i < element.num; i++) {
                    if (element.region_ID(i, 0) == k) {
                        IDs_to_request[count] = element.region_ID(i, 1);
                        IDs[count++] = i;

                        outstring << "    PROCESS " << proc_id << " WILL RECEIVE " << element.name
                                << " WITH LOCAL ID " << i << " FROM PROCESS " << k
                                << " (" << element.region_ID(i, 1) << ") " << std::endl;
                    }
                }

                if (count != num_to_receive) {
                    throw std::runtime_error(
                        "FindElementConnectivity: " + std::string(element.name) +
                        " count mismatch on process " + std::to_string(proc_id) +
                        " for region " + std::to_string(k));
                }

                element.receive_IDs[k] = IDs;

                MPI_Sendrecv(
                    IDs_to_request.data(), num_to_receive,  MPI_UNSIGNED, k, element.tag_ids,
                    IDs_to_send.data(),    num_to_send,     MPI_UNSIGNED, k, element.tag_ids,
                    MPI_COMM_WORLD, &status
                );

                for (unsigned i = 0; i < num_to_send; i++) {
                    outstring << "    PROCESS " << proc_id << " WILL SEND " << element.name
                            << " WITH LOCAL ID " << IDs_to_send[i] << " TO PROCESS " << k
                            << " (" << i << ")   " << std::endl;
                }

                // Reorder send IDs ascending, for better memory access when packing buffers.
                std::vector<unsigned> IDs_ordered(IDs_to_send);
                std::vector<unsigned> index_map(IDs_to_send.size());
                // index_map[ordered_index] gives the index in send_IDs[k].
                // Therefore:
                //     send_IDs_ordered[k][i] == send_IDs[k][send_IDs_map[k][i]]

                std::vector<std::pair<unsigned, unsigned>> pairs(IDs_to_send.size());
                for (unsigned i = 0; i < pairs.size(); ++i) {
                    pairs[i].first  = i;
                    pairs[i].second = IDs_ordered[i];
                }

                std::sort(pairs.begin(), pairs.end(),
                        [](auto & left, auto & right) { return left.second < right.second; });

                for (unsigned i = 0; i < pairs.size(); ++i) {
                    index_map[i]   = pairs[i].first;
                    IDs_ordered[i] = pairs[i].second;
                }

                element.send_IDs[k]         = IDs_to_send;
                element.send_IDs_ordered[k] = IDs_ordered;
                element.send_IDs_map[k]     = index_map;
            }
        }
    }
#endif



// Function sets up the relevant send and receive lists for each processor
int Mesh::FindInterconnectivity(void) {
#ifdef _MPI
    std::ostringstream outstring;

    MPIElementConnectivity node_conn {
        "NODE", node_num_ng, node_num,
        node_region_ID,
        recv_node_regions, send_node_regions,
        receive_node_IDs, send_node_IDs, send_node_IDs_ordered, send_node_IDs_map,
        TAG_NODE_COUNT, TAG_NODE_IDS
    };
    FindMPIElementConnectivity(node_conn, PROC_ID, NPROC, outstring);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef _DEBUG
    globals->Output->Write(OUT_MESSAGE, &outstring);
#endif

    MPIElementConnectivity face_conn {
        "FACE", face_num_ng, face_num,
        face_region_ID,
        recv_face_regions, send_face_regions,
        receive_face_IDs, send_face_IDs, send_face_IDs_ordered, send_face_IDs_map,
        TAG_FACE_COUNT, TAG_FACE_IDS
    };
    FindMPIElementConnectivity(face_conn, PROC_ID, NPROC, outstring);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef _DEBUG
    globals->Output->Write(OUT_MESSAGE, &outstring);
#endif

    MPIElementConnectivity vertex_conn {
        "VERTEX", vertex_num_ng, vertex_num,
        vertex_region_ID,
        recv_vertex_regions, send_vertex_regions,
        receive_vertex_IDs, send_vertex_IDs, send_vertex_IDs_ordered, send_vertex_IDs_map,
        TAG_VERTEX_COUNT, TAG_VERTEX_IDS
    };
    FindMPIElementConnectivity(vertex_conn, PROC_ID, NPROC, outstring);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef _DEBUG
    globals->Output->Write(OUT_MESSAGE, &outstring);
#endif


    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate request array with worst case scenario (should probably be bigger than this?)
    // mpi_requests.reserve(4 * NPROC); 
    request_frames.emplace_back();   // base frame — never popped, lives for the Mesh's lifetime

#endif
    return 1;
}

bool Mesh::CheckGhostContiguityForElement(const char * name,
                                           std::vector<std::vector<unsigned>> & receive_IDs)
{
    bool ok = true;

    for (unsigned k = 0; k < receive_IDs.size(); ++k) {
        const std::vector<unsigned> & IDs = receive_IDs[k];
        if (IDs.empty()) continue;

        for (size_t i = 1; i < IDs.size(); ++i) {
            if (IDs[i] != IDs[i-1] + 1) {
                std::cerr << "[RANK " << PROC_ID << "] CONTIGUITY VIOLATION in "
                          << name << " ghosts from rank " << k
                          << ": IDs[" << i-1 << "]=" << IDs[i-1]
                          << " -> IDs[" << i << "]=" << IDs[i]
                          << " (gap or out-of-order; expected " << IDs[i-1]+1 << ")"
                          << std::endl;
                ok = false;
            }
        }
    }

    return ok;
}

bool Mesh::CheckGhostContiguity(void)
{
    bool all_ok = 0;
#ifdef _MPI
 
    bool node_ok   = CheckGhostContiguityForElement("NODE",   receive_node_IDs);
    bool face_ok    = CheckGhostContiguityForElement("FACE",   receive_face_IDs);
    bool vertex_ok  = CheckGhostContiguityForElement("VERTEX", receive_vertex_IDs);

    all_ok = node_ok && face_ok && vertex_ok;

   // Make sure every rank's stdout/stderr has actually flushed and been seen
    // before any rank proceeds — makes the failure output easier to read when
    // running with mpirun, rather than interleaved arbitrarily.
    std::cerr.flush();
    MPI_Barrier(MPI_COMM_WORLD);

    int local_ok = all_ok ? 1 : 0;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
    all_ok = (global_ok != 0);


    if (PROC_ID == 0) {
        std::cerr << (all_ok
            ? "[CheckGhostContiguity] PASSED on all ranks — contiguous-block assumption holds.\n"
            : "[CheckGhostContiguity] FAILED on at least one rank — see per-rank output above.\n");
    }
#endif

    return all_ok;
}


Mesh::ExchangeScope::ExchangeScope(Mesh * mesh) : mesh_(mesh) {
    mesh_->PushExchangeFrame();
    mesh_->exchange_depth++;
}

Mesh::ExchangeScope::~ExchangeScope() {
    mesh_->exchange_depth--;
    mesh_->PopExchangeFrame();
}

void Mesh::PushExchangeFrame(void) {
#ifdef _MPI
    request_frames.emplace_back();
#ifdef _DEBUG
    if (exchange_depth >= 100) {
        throw std::runtime_error(
            "Mesh::PushExchangeFrame: exchange_depth (" + std::to_string(exchange_depth) +
            ") has reached the tag-spacing limit of 100. Increase the spacing between "
            "TAG_BASE_NODE_DATA/TAG_BASE_FACE_DATA/TAG_BASE_VERTEX_DATA/TAG_BASE_NODE_XYZ_DATA "
            "or reduce exchange nesting depth to avoid tag collisions.");
    }
#endif
#endif
}

void Mesh::PopExchangeFrame(void) {
#ifdef _MPI
    if (!request_frames.back().requests.empty()) {
        throw std::runtime_error(
            "Mesh::PopExchangeFrame: exchange scope closed with unresolved in-flight requests — missing a WaitAllExchanges() call before scope exit.");
    }
    request_frames.pop_back();
#endif
}

void Mesh::IrecvElement(double * field_base, 
                        unsigned stride,                           
                        std::vector<unsigned> & recv_regions,
                        std::vector<std::vector<unsigned>> & receive_IDs,
                        int tag_base)
{
#ifdef _MPI
    int tag = tag_base + exchange_depth;
    auto & frame = request_frames.back();
    for (unsigned k = 0; k < NPROC; ++k) {
        if (recv_regions[k] == 0) continue;

        unsigned size  = recv_regions[k] * stride;
        unsigned start = receive_IDs[k][0] * stride;   // contiguous-ghost-block assumption, as before

        frame.requests.emplace_back();
        MPI_Irecv(field_base + start, size, MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &frame.requests.back());
    }
#endif
}

void Mesh::IsendElement(double * field_base,
                        unsigned stride,
                        std::vector<unsigned> & send_regions,
                        std::vector<std::vector<unsigned>> & send_IDs_ordered,
                        std::vector<std::vector<unsigned>> & send_IDs_map,
                        int tag_base)
{
#ifdef _MPI
    int tag = tag_base + exchange_depth;
    auto & frame = request_frames.back();
    for (unsigned k = 0; k < NPROC; ++k) {
        if (send_regions[k] == 0) continue;

        unsigned size = send_regions[k];
        auto buf = std::make_unique<std::vector<double>>(size * stride);
        for (unsigned i = 0; i < size; ++i) {
            unsigned src = send_IDs_ordered[k][i] * stride;
            unsigned dst = send_IDs_map[k][i] * stride;
            for (unsigned c = 0; c < stride; ++c)
                (*buf)[dst + c] = field_base[src + c];
        }

        frame.requests.emplace_back();
        MPI_Isend(buf->data(), size * stride, MPI_DOUBLE, k, tag, MPI_COMM_WORLD, &frame.requests.back());
        frame.owned_buffers.push_back(std::move(buf));
    }
#endif
}

void Mesh::IrecvNode(Array1D<double> & field) {
#ifdef _MPI
    IrecvElement(&field(0), 1, recv_node_regions, receive_node_IDs, TAG_BASE_NODE_DATA);
#endif
}

void Mesh::IsendFace(Array1D<double> & field) {
#ifdef _MPI
    IsendElement(&field(0), 1, send_face_regions, send_face_IDs_ordered, send_face_IDs_map,
                       TAG_BASE_FACE_DATA);
#endif
}

void Mesh::IrecvFace(Array1D<double> & field) {
#ifdef _MPI
    IrecvElement(&field(0), 1, recv_face_regions, receive_face_IDs, TAG_BASE_FACE_DATA);
#endif
}

void Mesh::IsendNode(Array1D<double> & field) {
#ifdef _MPI
    IsendElement(&field(0), 1, send_node_regions, send_node_IDs_ordered, send_node_IDs_map,
                       TAG_BASE_NODE_DATA);
#endif
}

void Mesh::IrecvVertex(Array1D<double> & field) {
#ifdef _MPI
    IrecvElement(&field(0), 1, recv_vertex_regions, receive_vertex_IDs, TAG_BASE_VERTEX_DATA);
#endif
}

void Mesh::IsendVertex(Array1D<double> & field) {
#ifdef _MPI
    IsendElement(&field(0), 1, send_vertex_regions, send_vertex_IDs_ordered, send_vertex_IDs_map,
                       TAG_BASE_VERTEX_DATA);
#endif
}

void Mesh::IrecvNodeXYZ(Array2D<double> & field) {
#ifdef _MPI
    IrecvElement(&field(0,0), 3, recv_node_regions, receive_node_IDs, TAG_BASE_NODE_XYZ_DATA);
#endif
}
void Mesh::IsendNodeXYZ(Array2D<double> & field) {
#ifdef _MPI
    IsendElement(&field(0,0), 3, send_node_regions, send_node_IDs_ordered, send_node_IDs_map,
                       TAG_BASE_NODE_XYZ_DATA);
#endif
}

void Mesh::WaitAllExchanges(void) {
#ifdef _MPI
    auto & frame = request_frames.back();
    if (!frame.requests.empty()) {
        MPI_Waitall((int)frame.requests.size(), frame.requests.data(), MPI_STATUSES_IGNORE);
        frame.requests.clear();
    }
    frame.owned_buffers.clear();   // safe to free now — MPI has finished reading them
#endif
}

