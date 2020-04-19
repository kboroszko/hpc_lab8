/*
 * A template for the 2019 MPI lab at the University of Warsaw.
 * Copyright (C) 2016, Konrad Iwanicki.
 * Refactoring 2019, Łukasz Rączkowski
 */

#include <cassert>
#include <mpi.h>
#include "graph-base.h"
#include "graph-utils.h"

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

    Graph* dataToSend = nullptr;

    if(myRank == 0){
        dataToSend = allocateGraphPart(
                numVertices,
                0,
                numVertices
        );
    }

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
        int recipientRank = 1;
        int partStart = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank);
        int partEnd = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank + 1) - 1;
        for (int i = 0; i < graph->numVertices; ++i) {
            if(i < partStart){
                initializeGraphRow(graph->data[i], i, graph->numVertices);
            } else {
                initializeGraphRow(dataToSend->data[i], i, dataToSend->numVertices);
                MPI_Request request;
                MPI_Isend(dataToSend->data[i],
                        dataToSend->numVertices,
                        MPI_INT,
                        recipientRank,
                        0,
                        MPI_COMM_WORLD,
                        &request);
            }
            if(i == partEnd){
                recipientRank++;
                partStart = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank);
                partEnd = getFirstGraphRowOfProcess(numVertices, numProcesses, recipientRank + 1) - 1;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    } else {
        int rows = end - start;
        for(int i = 0; i < rows; i++){
            //recieve data synchronously
            MPI_Recv(dataToSend->data[i], dataToSend->numVertices, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

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

    /* FIXME: implement */
    for(int node=0; node<numProcesses; node++){
        MPI_Request request;
        MPI_Ibarrier(MPI_COMM_WORLD, &request);
        if(myRank == node){
            for(int i=0; i<rows; i++){
                printGraphRow(graph->data[i], 0, graph->numVertices);
            }
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
    }
}

void destroyGraph(Graph* graph, int numProcesses, int myRank) {
    freeGraphPart(graph);
}
