/* Standard blocker, part of Pogo - a finite element package to 
simulate elastic wave propagation on the GPU
Copyright (C) 2013 Peter Huthwaite

If you find Pogo useful in your academic work, please cite the relevant papers;
information on our latest papers is available at <http://www.pogo-fea.com/>.

This file is part of the blocker, and part of Pogo.

Pogo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Pogo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Pogo.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <list>

#include "block.h"
#include "node.h"
#include <cmath>

#include <sstream>
#include <string>

#include "../common/err.h"

//#include "mex.h"

bool block::checkLinks(vector<block> &blocks) {
    //check whether or not a block is adjacent to the ones it is supposed to be
    //correct if necessary
    list <int> linked;
    getLinks(linked);
    bool prevLinked = false;
    bool nextLinked = false;
    for (list<int>::iterator p = linked.begin();
        p != linked.end(); p++) {
        if (*p == prevBlock) {
            prevLinked = true;
        }
        if (*p == nextBlock) {
            nextLinked = true;
        }
    }
    bool retVal = true;
    if (prevBlock != -1 && !prevLinked) {
        //not attached, so remove this link
        blocks[prevBlock].nextBlock = -1;
        prevBlock = -1;
        retVal = false;
    }
    if (nextBlock != -1 && !nextLinked) {
        //not attached, so remove this link
        blocks[nextBlock].prevBlock = -1;
        nextBlock = -1;
        retVal = false;
    }
    return retVal;
}

bool block::deactivateRelink(vector<block> &blocks) {
    //check that next block is adjacent to previous block
    list <int> nextLinked;
    if (nextBlock < 0) {
        isActive = false;
        if (prevBlock > 0)
            blocks[prevBlock].nextBlock = -1;
        prevBlock = -1;
        return false;
    }
    blocks[nextBlock].getLinks(nextLinked);
    bool blocksLinked = false;
    for (list<int>::iterator p = nextLinked.begin();
         p != nextLinked.end(); p++) {
//        if (blockNum == 380)
//            cout << *p << endl;
        if (*p == prevBlock) {
            blocksLinked = true;
            break;
        }
    }
    if (blocksLinked) {
        cout << "Removing " << blockNum << " and linking adjacent: " << prevBlock << " and " << nextBlock << endl;
        isActive = false;
        blocks[nextBlock].prevBlock = prevBlock;
        blocks[prevBlock].nextBlock = nextBlock;
        nextBlock = -1;
        prevBlock = -1;
        return true;
    } else {
        cout << "Previous and next blocks aren't adjacent. Can't relink." << endl;
        isActive = false;
        if (nextBlock > 0)
            blocks[nextBlock].prevBlock = -1;
        if (prevBlock > 0)
            blocks[prevBlock].nextBlock = -1;
        nextBlock = -1;
        prevBlock = -1;
        return false;
    }
}


void block::renumber(int newNumber) {
    blockNum = newNumber;

    //renumber nodes in nodeBlocks
    for(list<node>::iterator p = intNodes.begin();
        p != intNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = blockNum;
    }

    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = blockNum;
    }
}

int block::join(vector<block> &blocks) {
    int joinBlock = -1;

    int bothNext = -1;
    int bothPrev = -1;

    if (nextBlock != -1) {
        joinBlock = nextBlock;
        bothNext = blocks[joinBlock].nextBlock;
        bothPrev = prevBlock;
    } else {
        joinBlock = prevBlock;
        if (joinBlock > 0)
            bothPrev = blocks[joinBlock].prevBlock;
        bothNext = nextBlock;
    }


    if (joinBlock == -1) {
        list <int> linkedBlocks;
        getLinks(linkedBlocks);

        for (list <int> :: iterator p = linkedBlocks.begin();
             p!= linkedBlocks.end(); p++) {
            if (blocks[*p].isActive) {
                joinBlock = *p;
                break;
            }
        }
        if (joinBlock != -1) {
            if (blocks[joinBlock].prevBlock == -1) {
                // joinBlock's previous slot is free, therefore we can go before
                bothPrev = -1;
                bothNext = blocks[joinBlock].nextBlock;
            } else if (blocks[joinBlock].nextBlock == -1) {
                // joinBlock's noxt slot is free, therefore we can go after
                bothPrev = blocks[joinBlock].prevBlock;
                bothNext = -1;
            } else {
                //both of the 'slots' for the adjacent block are used. No space.
                cout << "No joining block found. Adjacent one is fully linked, so can't be used." << endl;
                isActive = false;
                return -1;
            }
        } else {
            cout << "No joining block found." << endl;
            isActive = false;
            return -1;
        }
    }


    if (bothPrev == blockNum || bothPrev == joinBlock)
        bothPrev = -1;
    if (bothNext == blockNum || bothNext == joinBlock)
        bothNext = -1;


    //cout << "Sizes before: " << size() << "," << blocks[joinBlock - 1].size() << endl;
    modifiedTimes++;
    blocks[joinBlock].combineFrom(*this);

    cout << "Joining " << blockNum
         << " to " << joinBlock
         << "; current size of block vector: " << blocks.size() << endl;


    //#1 - update links for joinBlock
    if (bothNext == joinBlock)
        bothNext = -1;
    if (bothPrev == joinBlock)
        bothPrev = -1;
    blocks[joinBlock].nextBlock = bothNext;
    blocks[joinBlock].prevBlock = bothPrev;

    if (bothNext >= 0)
        blocks[bothNext].prevBlock = joinBlock;
    if (bothPrev >= 0)
        blocks[bothPrev].nextBlock = joinBlock;


    int lastBlock = blocks.size()-1;
    if (blockNum != lastBlock) {
        //#2 - move last block's data to this block and renumber

        //move last block to this block number
        modifiedTimes = 0;
        combineFrom(blocks[lastBlock]);

        nextBlock = blocks[lastBlock].nextBlock;
        prevBlock = blocks[lastBlock].prevBlock;
        if (nextBlock >= 0)
            blocks[nextBlock].prevBlock = blockNum;
        if (prevBlock >= 0)
            blocks[prevBlock].nextBlock = blockNum;
    } //otherwise everything is where it needs to be - no need to renumber

    blocks.pop_back();
    return joinBlock;
}

void block::combineFrom(block &b) {

    if (b.modifiedTimes > modifiedTimes)
        modifiedTimes = b.modifiedTimes;

    //renumber all the nodes in nodeBlocks first
    for(list<node>::iterator p = b.intNodes.begin();
        p != b.intNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = blockNum;
    }

    for(list<node>::iterator p = b.boundNodes.begin();
        p != b.boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = blockNum;
    }

    //now splice together the arrays
    intNodes.splice(intNodes.begin(),b.intNodes);
    boundNodes.splice(boundNodes.begin(),b.boundNodes);

    //and correct the boundaries
    tidyBoundaries();
}
void block::tidyBoundaries() {
    //correct any non-boundary nodes
    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        bool isBoundary = false;
        int currNodeNum = p->nodeNum;
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];
            if (linkNode < 0) {
                stringstream errSs;
                errSs << "BIG ERROR. linkMat has an undefined node (value " << linkNode << ") in node " << currNodeNum
                        << " at position " << linkCnt << " out of " << nNodesLinked[currNodeNum] << endl;
                err newError("Blocker - tidying boundaries", 2401, errSs.str());
                throw newError;
            }
            if (nodeBlocks[linkNode] != blockNum) {
                isBoundary = true;
                break;
            }
        }

        if (isBoundary == false) {
            flags[currNodeNum] = 1;
            intNodes.push_back(*p);
            p = boundNodes.erase(p);
            if (boundNodes.size() > 0)
                p--;
            else
                break;
        }
    }
}

void block::redoLists() {
    //clear the lists and start again
    boundNodes.clear();
    intNodes.clear();
    for (int cnt = 0; cnt < nNodes; cnt++) {
        if (nodeBlocks[cnt] == blockNum) {
            node newNode;
            newNode.nodeNum = cnt;
            boundNodes.push_back(newNode);
            flags[cnt] = 2;
        }
    }
    tidyBoundaries();
}

void block::activeNodes(list<node> &activeNodesList) {

    //return a list of active nodes (i.e. nodes connected to undefined areas)
    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        bool hasActiveNodes = false;
        int currNodeNum = p->nodeNum;
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];
            if (flags[linkNode] == 0) {
                hasActiveNodes = true;
                break;
            }
        }

        if (hasActiveNodes == true) {
            activeNodesList.push_back(*p);
        }
    }
}


int block::advanceGreedy(vector<block> &blocks,
             int maxNodesAdvance//set to -1 if unused, -2 if want to ignore block size too
                   ) {

    if (!isActive)
        return 0;

    int nAdded = 0;

    //loop through boundary nodes for connecting block
    list<node>::iterator p = boundNodes.begin();
    int nBoundNodes = boundNodes.size();
    for (int boundNodeCnt = 0; boundNodeCnt < nBoundNodes; boundNodeCnt++) {//fixed iterations
        int currNodeNum = p->nodeNum;
        if (flags[currNodeNum] == 5) {
            p++;
            continue;
        }

        //loop through linked nodes
        int * linkLine = &linkMat[(currNodeNum)*nMaxLinkedNodes];

        bool isBoundary = false;
        bool anyAdded = false;
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkLine[linkCnt];

            if (flags[linkNode] == 0) {

                flags[linkNode] = 2;
                nodeBlocks[linkNode] = blockNum;

                node addNode;
                addNode.nodeNum = linkNode;
                boundNodes.push_back(addNode);
                nAdded++;
                anyAdded = true;
                if ((nAdded >= maxNodesAdvance && maxNodesAdvance > 0) || (size() >= maxSize && maxNodesAdvance != -2)) {
                    break;
                }
            } else {
                if (nodeBlocks[linkNode] != blockNum) {
                    isBoundary = true;
                }
            }
        }

        if (!anyAdded) {
            flags[currNodeNum] = 5; //no point in going back to this
        }
        if (!isBoundary) {
            flags[currNodeNum] = 1;
            intNodes.push_back(*p);
            p = boundNodes.erase(p);
            if (boundNodes.size() > 0)
                p--;
            else
                break;
        }

        if ((nAdded >= maxNodesAdvance && maxNodesAdvance > 0) || (size() >= maxSize && maxNodesAdvance != -2)) {
            break;
        }

        p++;
    }

    //make sure there are no internal nodes marked as boundaries
    tidyBoundaries();//should be able to remove this

    //if (nAdded == 0) isActive = 0;
    lastAdvance = nAdded;
    return nAdded;
}


//Do not use this - ambiguous nodes are unused in the code
int block::advance(list<node> &ambiguousNodes, bool isGreedy, vector<block> &blocks,
             int maxNodesAdvance//set to -1 if unused, -2 if want to ignore block size too
                   ) {

    if (!isActive)
        return 0;

    int nAdded = 0;

    //loop through boundary nodes for connecting block
    list<node>::iterator p = boundNodes.begin();
    int nBoundNodes = boundNodes.size();
    for (int boundNodeCnt = 0; boundNodeCnt < nBoundNodes; boundNodeCnt++) {//fixed iterations
        int currNodeNum = p->nodeNum;

        //loop through linked nodes
        int * linkLine = &linkMat[(currNodeNum)*nMaxLinkedNodes];
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkLine[linkCnt];

            if (flags[linkNode] == 0) {
                bool singleLink = true;

                int * linkLineIn = &linkMat[(linkNode)*nMaxLinkedNodes];
                if (!isGreedy) {
                    //loop through linked nodes to check if it could be claimed by two
                    for (int linkCntIn = 0; linkCntIn < nNodesLinked[linkNode]; linkCntIn++) {
                        //int linkNodeIn = linkMat[(linkNode-1)*nMaxLinkedNodes+linkCntIn];
                        int linkNodeIn = linkLineIn[linkCntIn];
                        int linkNodeBlock = nodeBlocks[linkNodeIn];
                        bool linkIsActive = false;
                        if (linkNodeBlock >= 0)
                            linkIsActive = blocks[linkNodeBlock].isActive;
                        if (flags[linkNodeIn] != 0
                                && linkNodeBlock != blockNum
                                && linkIsActive) {
                            singleLink = false;
                            break;
                        }
                    }
                }
                if (singleLink) { //can only be claimed by one block, so we take it
                    flags[linkNode] = 2;
                    nodeBlocks[linkNode] = blockNum;

                    node addNode;
                    addNode.nodeNum = linkNode;
                    boundNodes.push_back(addNode);
                    nAdded++;
                    if ((nAdded >= maxNodesAdvance && maxNodesAdvance > 0) || (size() >= maxSize && maxNodesAdvance != -2)) {
//                        if (size() >= maxSize && maxNodesAdvance != -2) {
//                            cout << "Stopped advancing " << blockNum << " because maxNodesAdvance == " << maxNodesAdvance << endl;
//                        }
                        break;
                    }
                } else {
                    flags[linkNode] = 4;
                    node addNode;
                    addNode.nodeNum = linkNode;
                    bool isOriginal = true;
                    for (list<node>::iterator q = ambiguousNodes.begin(); q != ambiguousNodes.end(); q++) {
                        if (q->nodeNum == linkNode) {
                            isOriginal = false;
                            break;
                        }
                    }
                    if (isOriginal)
                        ambiguousNodes.push_back(addNode);
                }
            } else {
            }
        }

        //check if current node is now a boundary or not
        bool isBoundary = false;
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];

            if (nodeBlocks[linkNode] != blockNum) {
                isBoundary = true;
                break;
            }
        }
        if (!isBoundary) {
            flags[currNodeNum] = 1;
            intNodes.push_back(*p);
            p = boundNodes.erase(p);
            if (boundNodes.size() > 0)
                p--;
            else
                break;
        }

        if ((nAdded >= maxNodesAdvance && maxNodesAdvance > 0) || (size() >= maxSize && maxNodesAdvance != -2)) {
            break;
        }

        p++;
    }

    //make sure there are no internal nodes marked as boundaries
    tidyBoundaries();

    //if (nAdded == 0) isActive = 0;
    lastAdvance = nAdded;
    return nAdded;
}


void block::genNewBlock(vector<block> &blocks, block &newBlock) {
    //use this block to seed a new block

    if (!isActive)
        return;

    //block newBlock;
    newBlock.clear();
    newBlock.setup(*this);
    newBlock.blockNum = blocks.size();

    //loop through boundary nodes for connecting block
    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;

        //loop through linked nodes
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];

            if (flags[linkNode] == 0) {

                bool singleLink = true;

                //loop through linked nodes to that link to check if it could be claimed by two blocks
                for (int linkCntIn = 0; linkCntIn < nNodesLinked[currNodeNum]; linkCntIn++) {
                    int linkNodeIn = linkMat[(linkNode)*nMaxLinkedNodes+linkCntIn];
                    if (linkNodeIn >= 0
                            && flags[linkNodeIn] != 0
                            && nodeBlocks[linkNodeIn] != blockNum
                            && nodeBlocks[linkNodeIn] != newBlock.blockNum) {
                        singleLink = false;
                        //cout << nodeBlocks[linkNodeIn-1] << " vs " << connBlock+1 << ", " << blockCnt+1 << endl;
                        break;
                    }
                }
                if (singleLink) { //can only be claimed by one block, so we take it
                    flags[linkNode] = 2;
                    nodeBlocks[linkNode] = newBlock.blockNum;

                    node addNode;
                    addNode.nodeNum = linkNode;
                    newBlock.boundNodes.push_back(addNode);
                }
            }
        }
    }

    newBlock.parentBlock = blockNum;
    childBlock = newBlock.blockNum;
    if (newBlock.size() == 0) {
        childBlock = -1;
    }

    isActive = false;
}


void block::minimise(vector<block> &blocks) {
    //take the block back to its roots

    if (!isActive)
        return;

    //lose the internal nodes first
    for(list<node>::iterator p = intNodes.begin();
        p != intNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = -1;
        flags[currNodeNum] = 0;
    }
    intNodes.clear();

    //cout << "Done internals." << endl;
    //remove all boundary nodes except the ones linked to inactive blocks
    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;

        bool isInactiveLinked = false;

        //loop through linked nodes
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];
            int blockLink = nodeBlocks[linkNode];
            //int blockLink = blocks[linkNode - 1].blockNum;
            if (blockLink >= 0 && !blocks[blockLink].isActive)
                isInactiveLinked = true;
        }

        if (!isInactiveLinked) {
            nodeBlocks[currNodeNum] = -1;
            flags[currNodeNum] = 0;
            p = boundNodes.erase(p);
            if (boundNodes.size() > 0)
                p--;
            else
                break;
        } else {
            flags[currNodeNum] = 2;
        }
    }

}
void block::erase() {
    //remove all the data associated with a block
    for(list<node>::iterator p = intNodes.begin();
        p != intNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = -1;
        flags[currNodeNum] = 0;
    }
    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        nodeBlocks[currNodeNum] = -1;
        flags[currNodeNum-1] = 0;
    }
    boundNodes.clear();
    intNodes.clear();
}

bool block::removeFalseLinks(vector <block> &blocks) {
    list <int> linkedBlocks;
    getLinks(linkedBlocks);
    bool nextIsLinked = false;
    bool prevIsLinked = false;
    if (prevBlock == -1) prevIsLinked = true;
    if (nextBlock == -1) nextIsLinked = true;
    for (list<int>::iterator p = linkedBlocks.begin();
         p != linkedBlocks.end(); p++) {
        if (*p == nextBlock)
            nextIsLinked = true;
        if (*p == prevBlock)
            prevIsLinked = true;
    }
    if (!nextIsLinked) {
        blocks[nextBlock].prevBlock = -1;
        nextBlock = -1;
    }
    if (!prevIsLinked) {
        blocks[prevBlock].nextBlock = -1;
        prevBlock = -1;
    }
    if (!prevIsLinked || !nextIsLinked)
        return false;
    return true;
}

bool block::verify(vector<block> &blocks) {
    for(list<node>::iterator p = intNodes.begin();
        p != intNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        if (nodeBlocks[currNodeNum] != blockNum || flags[currNodeNum] == 0) {
            cout << "Internal nodes list not matching nodeBlocks array" << endl;
            return false;
        }
    }
    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;
        if (nodeBlocks[currNodeNum] != blockNum || flags[currNodeNum] == 0) {
            cout << "Boundary nodes list not matching nodeBlocks array" << endl;
            return false;
        }
        for(list<node>::iterator q = intNodes.begin();
            q != intNodes.end(); q++ ) {
            int intNodeNum = q->nodeNum;
            if (intNodeNum == currNodeNum) {
                cout << "A single node appears on both internal and boundary lists." << endl;
                cout << intNodeNum << endl;
                return false;
            }
        }
    }
    int totNodes = 0;
    for (int cnt = 0; cnt < nNodes; cnt++)
        if (nodeBlocks[cnt] == blockNum)
            totNodes++;

    if (totNodes != size()) {
        cout << "Total counts of lists and nodeBlocks array not matching" << endl;
        cout << "Block size: " << size() << " vs " << totNodes << " counted." << endl;
        return false;
    }

    if ((nextBlock < 0 && nextBlock != -1) ||
            (prevBlock < 0 && prevBlock != -1) ||
            abs(nextBlock) > blocks.size()-1 ||
            abs(prevBlock) > blocks.size()-1) {
        cout << "nextBlock and prevBlock not in range: " << nextBlock << "," << prevBlock << endl;
        cout << "Size: " << blocks.size() << endl;
        return false;
    }
    if (nextBlock != -1 && blockNum != blocks[nextBlock].prevBlock) {
        cout << "prevBlock (" << blocks[nextBlock].prevBlock
             << ") of nextBlock ("
             << nextBlock << ") doesn't point to us..." << endl;
        return false;
    }
    if (prevBlock != -1 && blockNum != blocks[prevBlock].nextBlock) {
        cout << "nextBlock (" << blocks[prevBlock].nextBlock
             << ") of prevBlock (" << prevBlock << ") doesn't point to us..." << endl;
        return false;
    }

    list <int> linkedBlocks;
    getLinks(linkedBlocks);
    bool nextIsLinked = false;
    bool prevIsLinked = false;
    if (prevBlock == -1) prevIsLinked = true;
    if (nextBlock == -1) nextIsLinked = true;
    for (list<int>::iterator p = linkedBlocks.begin();
         p != linkedBlocks.end(); p++) {
        if (*p == nextBlock)
            nextIsLinked = true;
        if (*p == prevBlock)
            prevIsLinked = true;
    }
    if (!nextIsLinked)
        cout << "nextBlock " << nextBlock << " is not physically linked." << endl;
    if (!prevIsLinked)
        cout << "prevBlock " << prevBlock << " is not physically linked." << endl;
    if (!prevIsLinked || !nextIsLinked)
        return false;


    return true;
}

//void block::fillToDepth(double * nodeDepth, int dLim) {
//    int nNodesInBlock = size();

//    for (int cnt = 0; cnt < nNodes*2; cnt++) {
//        node currNode = boundNodes.front();
//        int currNodeNum = currNode.nodeNum;


//        bool isDepthBoundary = false;
//        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
//            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];
//            if (nodeDepth[linkNode] > dLim) {
//                isDepthBoundary = true;
//            } else {
//                if (flags[linkNode] == 0) {
//                    node addNode;
//                    addNode.nodeNum = linkNode;
//                    addNode.score = currNode.score-1;
//                    boundNodes.push_back(addNode);

//                    flags[linkNode] = 2;
//                    nodeBlocks[linkNode] = blockNum;
//                    nNodesInBlock++;
//                }
//            }
//        }


//        if (isDepthBoundary == false) {
//            intNodes.push_back(currNode);
//            boundNodes.pop_front();
//            flags[currNode.nodeNum] = 1;
//        } else {
//            //put the current node to the end
//            currNode.score = -1;
//            boundNodes.push_back(currNode);
//            boundNodes.pop_front();
//            flags[currNode.nodeNum] = 3;
//        }
//        if (nNodesInBlock >= maxSize) break;
//    }

////    tidy up
//    tidyBoundaries();
//    //block.isActive = false;

//}
void block::getLinks(list<int> &linkedBlocks) {

    linkedBlocks.clear();

    for(list<node>::iterator p = boundNodes.begin();
        p != boundNodes.end(); p++ ) {
        int currNodeNum = p->nodeNum;

        //loop through linked nodes
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
            int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];
            int blockLink = nodeBlocks[linkNode];
            if (blockLink >= 0 && blockLink != blockNum) {
                //check to see if it's already in the list
                bool alreadyExists = false;
                for(list<int>::iterator q = linkedBlocks.begin();
                    q != linkedBlocks.end(); q++ ) {
                    if (*q == blockLink) {
                        alreadyExists = true;
                        break;
                    }
                }
                if (!alreadyExists) {
                    linkedBlocks.push_back(blockLink);
                }
            }
        }
    }
}



void block::depthInBlock(int * depths, int &deepVal, int &deepNode) {
    int currDepth = 1;
    deepVal = 1;
    deepNode = -1;

    for (int cnt = 0; cnt < nNodes*2; cnt++) { //just use a for loop to avoid possible infinite loop
        int nextDepth = currDepth + 1;
        int updatePoints = 0;
        //for (int nCnt = 0; nCnt < nNodes; nCnt++) {//more efficient just to use list of nodes in block
        for (list <node>::iterator p = intNodes.begin(); p != intNodes.end(); p++) {
            int nCnt = p->nodeNum;
            if (depths[nCnt] == currDepth) {
                int nNLinkHere = nNodesLinked[nCnt];
                for (int connCnt = 0; connCnt < nNLinkHere; connCnt++) {
                    int checkNode = linkMat[nCnt*nMaxLinkedNodes+connCnt];
                    if (depths[checkNode] == 0 && nodeBlocks[checkNode] == blockNum) {
                        depths[checkNode] = nextDepth;
                        deepNode = checkNode;
                        updatePoints++;
                    }
                }
            }
        }
        for (list <node>::iterator p = boundNodes.begin(); p != boundNodes.end(); p++) {
            int nCnt = p->nodeNum;

            if (depths[nCnt] == currDepth) {
                int nNLinkHere = nNodesLinked[nCnt];
                for (int connCnt = 0; connCnt < nNLinkHere; connCnt++) {
                    int checkNode = linkMat[nCnt*nMaxLinkedNodes+connCnt];
                    if (depths[checkNode] == 0 && nodeBlocks[checkNode] == blockNum) {
                        depths[checkNode] = nextDepth;
                        deepNode = checkNode;
                        updatePoints++;
                    }
                }
            }
        }

        if (updatePoints == 0) break;
        currDepth = nextDepth;
    }
    deepVal = currDepth;
}

void block::separateUnjoined(vector<block> &blocks) { //need to sort out the stuff with the linked blocks
    //does a flood fill to check whether a block is joined. If not, it separates them.
    int * depths = new int[nNodes];

    int localBlockNum = blockNum;
    block * thisBlock = this;
    // this remains unchanged. When we start using push_back on blocks
    //(which contains this block) it can rearrange the memory locations for all the data - we
    //therefore need to refer explicitly to any contained data.

    int currNode = boundNodes.front().nodeNum;

    for (int cnt = 0; cnt < nNodes; cnt++)
        depths[cnt] = 0;

    depths[currNode] = 1;
    int maxDepth, maxNode;
    depthInBlock(depths, maxDepth, maxNode);

    int nDeep = 0;
    for (int cnt = 0; cnt < nNodes; cnt++)
        if (depths[cnt] > 0)
            nDeep++;

    bool blockChanged = false;

    //now loop through the block
    //for (int cnt = 0; cnt < thisBlock->nNodes; cnt++) {
    for (list<node>::iterator p = boundNodes.begin(); p != boundNodes.end(); p++) {
        int cnt = p->nodeNum;
        if (depths[cnt] == 0) {
            blockChanged = true;
            //found one in a separate block
            //sepSectionToNewBlock(thisBlock, blocks, cnt, localBlockNum, depths );
        //if (thisBlock->nodeBlocks[cnt] == thisBlock->blockNum && depths[cnt] == 0) {
//            //cout << "Separated block found: " << blockNum << endl;//mexEvalString("drawnow;");

            block newBlock;
            newBlock.setup(*thisBlock);
            //newBlock.setup(nNodes, nMaxLinkedNodes, nNodesLinked, linkMat, flags, nodeBlocks);
            newBlock.blockNum = blocks.size();

            //seed it
            node addNode;
            addNode.nodeNum = cnt;
            newBlock.boundNodes.push_back(addNode);
            thisBlock->nodeBlocks[cnt] = newBlock.blockNum;
            thisBlock->flags[cnt] = 2;

            while (1) {
                list<node>::iterator p = newBlock.boundNodes.begin();
                int nBoundNodes = newBlock.boundNodes.size();
                int nAdded = 0;
                for (int boundNodeCnt = 0; boundNodeCnt < nBoundNodes; boundNodeCnt++) {
                    int currNodeNum = p->nodeNum;

                    for (int linkCnt = 0; linkCnt < thisBlock->nNodesLinked[currNodeNum]; linkCnt++) {//loop through linked nodes
                        int linkNode = thisBlock->linkMat[(currNodeNum)*thisBlock->nMaxLinkedNodes+linkCnt];

                        if (thisBlock->nodeBlocks[linkNode] == thisBlock->blockNum) {
                            node addNode;
                            addNode.nodeNum = linkNode;
                            newBlock.boundNodes.push_back(addNode);
                            thisBlock->nodeBlocks[linkNode] = newBlock.blockNum;
                            thisBlock->flags[linkNode] = 2;
                            depths[linkNode] = 1;
                            nAdded++;
                        }
                    }
                    p++;
                }
                //cout << "nAdded: " << nAdded << endl;mexEvalString("drawnow;");
                if (nAdded == 0)
                    break;
            }

            newBlock.tidyBoundaries();
            newBlock.prevBlock = -1;
            newBlock.nextBlock = -1;
            list <int> linkedBlocks;
            newBlock.getLinks(linkedBlocks);
            for (list<int>::iterator p = linkedBlocks.begin();
                 p != linkedBlocks.end(); p++) {
                if (*p == thisBlock->prevBlock) {
                    newBlock.prevBlock = *p;
                    blocks[prevBlock].nextBlock = newBlock.blockNum;
                    thisBlock->prevBlock = -1;
                }
                if (*p == thisBlock->nextBlock) {
                    newBlock.nextBlock = *p;
                    blocks[nextBlock].prevBlock = newBlock.blockNum;
                    thisBlock->nextBlock = -1;
                }
            }

            cout << "Splitting block - new block number: " << newBlock.blockNum
                 << "(original is " << localBlockNum << ")" << endl;
            blocks.push_back(newBlock);
            thisBlock = &(blocks[localBlockNum]);
        }
    }

    delete [] depths;

    if (blockChanged)
        thisBlock->redoLists();
}


void block::findFurthest(int deepLoc[2], int &currNode) {
    deepLoc[0] = -1;
    deepLoc[1] = -1;
    bool isGood[2];
    int * depths = new int[nNodes];

    for (int cnt = 0; cnt < 11; cnt++) {
        for (int nCnt = 0; nCnt < nNodes; nCnt++) {
            depths[nCnt] = 0;
        }
        depths[currNode] = 1;
        int maxDepth, deepestNode;
        depthInBlock(depths, maxDepth, deepestNode);
        int maxLoc = -1;
        maxLoc = deepestNode;
        if (maxLoc == -1) {
            stringstream errSs;
            errSs << "No maximum depth found in block " << blockNum << endl;
            err newError("Blocker - block splitting", 2403, errSs.str());
            throw newError;
        }
        if (maxDepth == 1) {
            stringstream errSs;
            errSs << "Maximum depth == 1, the starting depth. Something has gone wrong." << endl;
            errSs << "In block " << blockNum << endl;
            err newError("Blocker - block splitting", 2404, errSs.str());
            throw newError;
        }
        if (deepLoc[cnt%2] == maxLoc) {
            isGood[cnt%2] = true;
            if (isGood[(cnt+1)%2]) { //both nodes are being repeated
                cout << cnt << endl;
                break;
            }
        } else
            isGood[cnt%2] = false;
        deepLoc[cnt%2] = maxLoc;
        currNode = maxLoc;
    }
    delete [] depths;
}

int block::split(vector<block> &blocks) {
    //splits a block in 2, starting from deepest points
    //returns number of new block

    if (size() == 1) {
        cout << "Can't split block - it's only got one node." << endl;
        return -1;
    }


    //do depth projection first
    int currNode = -1;
    for (list <node> :: iterator p = boundNodes.begin();
         p!= boundNodes.end(); p++) {
        currNode = p->nodeNum;
//        if (blockNum == 2139) {
//            cout << "Testing " << currNode << endl;
//        }
        bool hasLinks = false;
        for (int linkCnt = 0; linkCnt < nNodesLinked[currNode]; linkCnt++) {//loop through linked nodes
            int linkNode = linkMat[(currNode)*nMaxLinkedNodes+linkCnt];
            if (nodeBlocks[linkNode] == blockNum) {
                hasLinks = true;
                break;
            }
        }
        if (hasLinks)
            break;
    }
    if (currNode == -1) {
        stringstream errSs;
        errSs << "Unable to find a good starting node in block " << blockNum << endl;
        err newError("Blocker - block splitting", 2402, errSs.str());
        throw newError;
    }

    int deepLoc[2];
    //cout << "finding deepest" << endl;
    findFurthest(deepLoc, currNode);
    //cout << "done" << endl;

    //now deepLoc[] contains the two 'furthest' nodes

    //cout << "Deepest: " << deepLoc[0]+1 << "," << deepLoc[1]+1 << endl;
    if (deepLoc[0] == deepLoc[1]) {
        stringstream errSs;
        errSs << "Deepest points are the same, in " << blockNum << endl;
        err newError("Blocker - block splitting", 2405, errSs.str());
        throw newError;
    }

    //define new blocks
    modifiedTimes++;

    block nb[2];
    nb[1].blockNum = blocks.size() + 0;
    nb[0].blockNum = blocks.size() + 1;
    nb[0].setup(*this);
    nb[1].setup(*this);
    nb[0].modifiedTimes = modifiedTimes;
    nb[1].modifiedTimes = modifiedTimes;



    for (int cnt = 0; cnt < 2; cnt++) {
        node addNode;
        addNode.nodeNum = deepLoc[cnt];
        nb[cnt].boundNodes.push_back(addNode);
        nodeBlocks[deepLoc[cnt]] = nb[cnt].blockNum;
        flags[deepLoc[cnt]] = 2;
    }

    int nItsMax = size(); //upper bound so we don't get infinite loop



    bool tryOther = false;
    for (int cnt = 0; cnt < nItsMax; cnt++) {//loop through rows
        int nAdded = 0;
        int sCnt;
        if ((nb[0].size() > nb[1].size()) != tryOther)
            sCnt = 1;
        else
            sCnt = 0;


            list<node>::iterator p = nb[sCnt].boundNodes.begin();
            int nBoundNodes = nb[sCnt].boundNodes.size();
            for (int boundNodeCnt = 0; boundNodeCnt < nBoundNodes; boundNodeCnt++) {
                int currNodeNum = p->nodeNum;

                for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {//loop through linked nodes
                    int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];

                    if (nodeBlocks[linkNode] == blockNum) {
                        node addNode;
                        addNode.nodeNum = linkNode;
                        nb[sCnt].boundNodes.push_back(addNode);
                        nodeBlocks[linkNode] = nb[sCnt].blockNum;
                        flags[linkNode] = 2;
                        nAdded++;
                    }
                }
                p++;
            }
        if (nAdded == 0) {
            if (tryOther)
                break;
            else
                tryOther = true;
        }
    }


    nb[0].redoLists();
    nb[1].redoLists();





    int beforeAfter[2]; // = 1 for before, = 2 for after, = 0 for undefined

    list<int> adjacentBlocks;
    for (int cnt = 0; cnt < 2; cnt++) {
        beforeAfter[cnt] = 0;

        nb[cnt].getLinks(adjacentBlocks);

        for (list<int>::iterator p = adjacentBlocks.begin();
             p != adjacentBlocks.end(); p++) {
            if (*p == prevBlock && beforeAfter[cnt] % 2 == 0) {
                beforeAfter[cnt] += 1;
            }
            if (*p == nextBlock && ((beforeAfter[cnt]/2) % 2) == 0) {
                beforeAfter[cnt] += 2;
            }
        }
        adjacentBlocks.clear();
    }



    int ord = 0; //= 1 put 0 first, = 2 put 1 first.

    //decide order based on links. 0 dominates if they both want the same.
    if (beforeAfter[0] == 0 && beforeAfter[1] == 0) {
        cout << "Can't find adjacent blocks " << blockNum << endl;
    } else if (beforeAfter[0] == 3 || beforeAfter[0] == 0){ // 0 no preference
        if (beforeAfter[1] == 3) {
            //either way - it doesn't matter
            ord = 1;
        } else if (beforeAfter[1] == 2) { //1 prefers to go after
            ord = 1;
        } else if (beforeAfter[1] == 1) { //1 prefers to go first
            ord = 2;
        } else { //1 doesn't link. Order doesn't matter.
            ord = 1;
        }
    } else if (beforeAfter[0] == 2) { // 0 prefers to go second
        ord = 2;
    } else if (beforeAfter[0] == 1) { // 0 prefers to go first
        ord = 1;
    }


    //set one block to replace the original
    nb[0].renumber(blockNum);


    //we may have missed an unlinkable node/section
    redoLists();

    int bothPrev = prevBlock;
    int bothNext = nextBlock;

    if (ord == 1) {
        //this first
        prevBlock = bothPrev;
        nextBlock = nb[1].blockNum;

        //second
        nb[1].prevBlock = blockNum;
        nb[1].nextBlock = bothNext;

        if (bothPrev >= 0)
            blocks[bothPrev].nextBlock = blockNum;
        if (bothNext >= 0)
            blocks[bothNext].prevBlock = nb[1].blockNum;

    } else {
        //this second
        prevBlock = nb[1].blockNum;
        nextBlock = bothNext;

        //first
        nb[1].prevBlock = bothPrev;
        nb[1].nextBlock = blockNum;

        if (bothPrev >= 0)
            blocks[bothPrev].nextBlock = nb[1].blockNum;
        if (bothNext >= 0)
            blocks[bothNext].prevBlock = blockNum;
    }


    cout << "New block num: " << nb[1].blockNum << ",";
    cout << "New block size: " << nb[1].size() << "," << size() << endl;

    //set the other to a new one at the end
    blocks.push_back(nb[1]);


    return blocks.size()-1;//return the number of the new block
}

int block::splitInLayer(vector<block> &blocks, bool limitSize) {
    //splits a long block into 2

    if (limitSize && size() > XBLOCKSIZE*YBLOCKSIZE*2) {
        cout << "Block is too large to be split into two." << endl;
        cout << "Size: " << size() << ", block: " << blockNum << endl;
        return -1;
    }

    if (prevBlock == -1 && nextBlock == -1) {
        //not attached on either side, so do standard split.
        return split(blocks);
    }

    modifiedTimes++;

    //define new blocks
    block nb[2];
    nb[1].blockNum = blocks.size() + 0;
    nb[0].blockNum = blocks.size() + 1;
    nb[0].setup(*this);
    nb[1].setup(*this);
    nb[0].modifiedTimes = modifiedTimes;
    nb[1].modifiedTimes = modifiedTimes;

    int boundBlocks[2];
    boundBlocks[0] = prevBlock;
    boundBlocks[1] = nextBlock;



    //seed new blocks
    for (int cnt = 0; cnt < 2; cnt++) {
        if (boundBlocks[cnt] != -1) {
            for (list<node>::iterator p = blocks[boundBlocks[cnt]].boundNodes.begin();
                 p != blocks[boundBlocks[cnt]].boundNodes.end(); p++ ) {
                int currNodeNum = p->nodeNum;
                for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {
                    int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];
                    int linkBlockNum = nodeBlocks[linkNode];
                    if (linkBlockNum == blockNum) {
                        node addNode;
                        addNode.nodeNum = linkNode;
                        nb[cnt].boundNodes.push_back(addNode);
                        nodeBlocks[linkNode] = nb[cnt].blockNum;
                        flags[linkNode] = 2;
                    }
                }
            }
        }
    }
    if (nb[0].size() == 0 && nb[1].size() == 0) {
        stringstream errSs;
        errSs << "Both blocks seeded with zero size. Block " << blockNum << endl;
        errSs << nextBlock << "," << prevBlock << endl;
        err newError("Blocker - block splitting in layer", 2406, errSs.str());
        throw newError;
    }

    for (int cnt = 0; cnt < 2; cnt++) {
        if (nb[cnt].size() == 0) {
            //in this case we seed using the deepest point from the existing adjacent block
            int * depths = new int[nNodes];
            for (int nCnt = 0; nCnt < nNodes; nCnt++)
                depths[nCnt] = 0;

            for (list<node>::iterator p = nb[(cnt+1)%2].boundNodes.begin();
                 p != nb[(cnt+1)%2].boundNodes.end(); p++) {
                depths[(*p).nodeNum] = 1;
            }
            int maxDepth, maxLoc;
            depthInBlock(depths, maxDepth, maxLoc);

            if (maxLoc == -1) {
                stringstream errSs;
                errSs << "No maximum depth found. Block " << blockNum << endl;
                err newError("Blocker - block splitting in layer", 2407, errSs.str());
                throw newError;
            }

            node addNode;
            addNode.nodeNum = maxLoc;
            nb[cnt].boundNodes.push_back(addNode);
            nodeBlocks[maxLoc] = nb[cnt].blockNum;
            flags[maxLoc] = 2;

            delete [] depths;
        }
    }

    if (nb[0].size() == 0 || nb[1].size() == 0) {
        stringstream errSs;
        errSs << "One block seeded with zero size. Block " << blockNum << endl;
        errSs << nextBlock << "," << prevBlock << endl;
        err newError("Blocker - block splitting in layer", 2408, errSs.str());
        throw newError;
    }


    int nItsMax = size(); //upper bound so we don't get infinite loop

    //advance each, doing smallest block first, keeping within original block
    bool tryOther = false;
    for (int cnt = 0; cnt < nItsMax; cnt++) {//loop through rows
        int nAdded = 0;
        int sCnt;
        if ((nb[0].size() > nb[1].size()) != tryOther)
            sCnt = 1;
        else
            sCnt = 0;
        //for (int sCnt = 0; sCnt < 2; sCnt++) {// do each block
            if (limitSize && nb[sCnt].size() >= nb[sCnt].maxSize)
                continue;

            list<node>::iterator p = nb[sCnt].boundNodes.begin();
            int nBoundNodes = nb[sCnt].boundNodes.size();
            for (int boundNodeCnt = 0; boundNodeCnt < nBoundNodes; boundNodeCnt++) {
                int currNodeNum = p->nodeNum;

                for (int linkCnt = 0; linkCnt < nNodesLinked[currNodeNum]; linkCnt++) {//loop through linked nodes
                    int linkNode = linkMat[(currNodeNum)*nMaxLinkedNodes+linkCnt];

                    if (nodeBlocks[linkNode] == blockNum) {
                        node addNode;
                        addNode.nodeNum = linkNode;
                        nb[sCnt].boundNodes.push_back(addNode);
                        nodeBlocks[linkNode] = nb[sCnt].blockNum;
                        flags[linkNode] = 2;
                        nAdded++;
                    }
                }
                p++;
            }
        if (nAdded == 0) {
            if (tryOther)
                break;
            else
                tryOther = true;
        }
    }

    nb[0].tidyBoundaries();
    nb[1].tidyBoundaries();

    //set one block to replace the original
    nb[0].renumber(blockNum);
//    boundNodes.clear();
//    intNodes.clear();

//    intNodes.splice(intNodes.begin(),nb[0].intNodes);
//    boundNodes.splice(boundNodes.begin(),nb[0].boundNodes);
    redoLists();

    //do links
    int bothPrev = prevBlock;
    int bothNext = nextBlock;
    prevBlock = bothPrev;
    nextBlock = nb[1].blockNum;
    nb[1].prevBlock = blockNum;
    nb[1].nextBlock = bothNext;
    if (bothPrev >= 0)
        blocks[bothPrev].nextBlock = blockNum;
    if (bothNext >= 0)
        blocks[bothNext].prevBlock = nb[1].blockNum;

    cout << "New block num: " << nb[1].blockNum << ",";
    cout << "New block size: " << nb[1].size() << ", current block size: " << size() << endl;


    //set the other to a new one at the end
    blocks.push_back(nb[1]);

    return blocks.size()-1;//return the number of the new block
}

int block::getStatus(vector<block> &blocks) {
    if (isActive) {
        blockFlag = 0;
    } else {
        blockFlag = 2; // i.e. all linked blocks are inactive (check in a second)
        list <int> linkedBlocks;
        getLinks(linkedBlocks);
        for (list <int>::iterator p = linkedBlocks.begin();
             p != linkedBlocks.end(); p++) {
            if (blocks[*p].isActive) {
                blockFlag = 1;
                break;
            }
        }
    }
    return blockFlag;
}


void verifyAll(vector<block> &blocks) {
    for (unsigned int cnt = 0; cnt < blocks.size(); cnt++)
        blocks[cnt].verify(blocks);
}

void separate(int nNodes, int nMaxLinkedNodes, int *nNodesLinked, int *linkMat, short *flags, int *nodeBlocks, vector<block> &blocks) {
    //using nodeBlocks[], populate the blocks vector so it can be used for memory linking
    blocks.clear();

    int maxBlock = 0;
    for (int cnt = 0; cnt < nNodes; cnt++) {
        if (maxBlock < nodeBlocks[cnt])
            maxBlock = nodeBlocks[cnt];
    }

    for (int currBlock = 0; currBlock < maxBlock+1; currBlock++) {
        block newBlock;
        newBlock.setup(nNodes, nMaxLinkedNodes, nNodesLinked, linkMat, flags, nodeBlocks);
        newBlock.blockNum = currBlock;
        blocks.push_back(newBlock);
    }

    int nBlocks = blocks.size();

    for (int cnt = 0; cnt < nNodes; cnt++) {
        if (nodeBlocks[cnt] >= 0) {
            if (nodeBlocks[cnt] > nBlocks-1) {
                stringstream errSs;
                errSs << "Something went wrong separating the block definition array into blocks" << endl;
                errSs << nodeBlocks[cnt] << " block number found vs number of blocks: " << nBlocks;
                err newError("Blocker - getting block definitions", 2409, errSs.str());
                throw newError;
            }
            node newNode;
            newNode.nodeNum = cnt;
            blocks[nodeBlocks[cnt]].boundNodes.push_back(newNode);
            flags[cnt] = 2;
        }
    }

    for (int currBlock = 0; currBlock < maxBlock+1; currBlock++) {
        blocks[currBlock].tidyBoundaries();
    }
}
