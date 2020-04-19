/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <cassert>
#include <mpi.h>
#include "graph-base.h"
#include "graph-utils.h"
#include <malloc.h>

int getFirstGraphRowOfProcess(int numVertices, int numProcesses, int myRank) {
    if(myRank == numProcesses){
        return numVertices;
    }
    int rowsInOne = (numVertices + numProcesses - 1)/numProcesses;
    return myRank * rowsInOne;
}

Graph* createAndDistributeGraph(int numVertices, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);

    Graph* graph = nullptr;
    int start = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank);
    int end = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1);


    graph = allocateGraphPart(
            numVertices,
            start,
            end
    );


    if (graph == nullptr) {
        return nullptr;
    }

    assert(graph->numVertices > 0 && graph->numVertices == numVertices);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    if(myRank == 0){
//        int buf_size = numVertices*rowsInOne*(numProcesses-1) + MPI_BSEND_OVERHEAD;
//        int * b = new int[buf_size];
//        MPI_Buffer_attach( b, buf_size*sizeof(int) );

        int recipientRank = 1;
        int partStart = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank);
        int partEnd = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank + 1) - 1;

        int *row = new int[numProcesses];

        for (int i = 0; i < graph->numVertices; ++i) {
            if(i < partStart){
                initializeGraphRow(graph->data[i], i, graph->numVertices);
            } else {
                initializeGraphRow(row, i, numVertices);
//                MPI_Request request;
                MPI_Send(row,
                        numVertices,
                        MPI_INT,
                        recipientRank,
                        0,
                        MPI_COMM_WORLD
//                        ,&request
                        );
            }
            if(i == partEnd){
                recipientRank++;
                partStart = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank);
                partEnd = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank + 1) - 1;
            }
        }

        delete[] row;
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        int rows = end - start;
        for(int i = 0; i < rows; i++){
            //recieve data synchronously
            MPI_Recv(graph->data[i], graph->numVertices, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

//    freeGraphPart(dataToSend);

    //  node 0:
    //      create whole graph
    //      send rows to everyone
    //      free whole graph
    //      wait for everyone to finish?
    //
    //  node other than 0:
    //  receive my part of the graph
    //  wait for everyone to finish?
    return graph;
}

void collectAndPrintGraph(Graph* graph, int numProcesses, int myRank) {
    assert(numProcesses >= 1 && myRank >= 0 && myRank < numProcesses);
    assert(graph->numVertices > 0);
    assert(graph->firstRowIdxIncl >= 0 && graph->lastRowIdxExcl <= graph->numVertices);

    int numVertices = graph->numVertices;
    int start = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank);
    int end = getFirstGraphRowOfProcess(numVertices, numProcesses, myRank + 1);
    int rows = end - start;

    int* recv_data = nullptr;
    int* send_data = nullptr;

    int rowsInOne = (numVertices + numProcesses - 1)/numProcesses;

    if(myRank == 0){
        recv_data = new int[(numVertices*rowsInOne*(numProcesses-1))];
    } else {
        send_data = new int[numVertices*rowsInOne];
        for(int i=0; i<rows; i++){
            for(int j=0; j < numVertices; j++){
                send_data[i*numVertices + j] = graph->data[i][j];
            }
        }
    }

    MPI_Gather(
            send_data,
            numVertices*rowsInOne,
            MPI_INT,
            recv_data,
            numVertices*rowsInOne,
            MPI_INT,
            0,
            MPI_COMM_WORLD);

    if(myRank == 0){
        for(int i=0; i<rows; i++)
            printGraphRow(graph->data[i], 0, graph->numVertices);
        for(int i=0; i<(graph->numVertices - rows); i++){
            printGraphRow(recv_data + (i*numVertices),0, numVertices );
        }
        delete[] recv_data;
    } else {
        delete[] send_data;
    }
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
