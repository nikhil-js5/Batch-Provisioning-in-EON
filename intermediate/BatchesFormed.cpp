#include <iostream>
#include <vector>
#include <cstdlib> // for rand()
#include <ctime>   // for srand()
#include <limits>  // for numeric_limits
#include <queue>   // for priority_queue
#include <algorithm> // for std::sort
#include <random>

using namespace std;

class LightpathRequest {
public:
    int getId() const {
        return id;
        }
    int getSourceNode() const {
        return sourceNode;
        }
    int getDestinationNode() const {
        return destinationNode;
        }
    int getSlotsRequired() const {
        return slotsRequired;
        }
    int getNumHops() const {
        return numHops;
        }
    bool isConnectionPossible() const {
        return connectionPossible;
        }

    void setId(int newId) {
        id = newId;
        }
    void setSourceNode(int newSourceNode) {
        sourceNode = newSourceNode;
        }
    void setDestinationNode(int newDestinationNode) {
        destinationNode = newDestinationNode;
        }
    void setSlotsRequired(int newSlotsRequired) {
        slotsRequired = newSlotsRequired;
        }
    void setNumHops(int newNumHops) {
        numHops = newNumHops;
        }
    void setConnectionPossible(bool newConnectionPossible) {
        connectionPossible = newConnectionPossible;
        }

private:
    int id;
    int sourceNode;
    int destinationNode;
    int slotsRequired;
    int numHops;

    bool connectionPossible;
};

//ALL FUNCTIONS IN THE CODE
vector<vector<int>> createNetwork(int numNodes);
void calculateShortestPath(const vector<vector<int>>& network, int sourceNode, vector<int>& dist);
int calculateTotalSlotsRequired(const std::vector<LightpathRequest>& batch);
std::vector<std::vector<LightpathRequest>> batchFormation(const std::vector<LightpathRequest>& superSet, int batchCapacity);


vector<vector<int>> createNetwork(int numNodes) {
    vector<vector<int>> network(numNodes, vector<int>(numNodes, 0));

    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            int isConnected;
            if ((std::rand() % 100 + 1) <= 40) {
                isConnected = 1;
            } else {
                isConnected = 0;
            }

            if (isConnected) {
                network[i][j] = 1;
                network[j][i] = 1;
            }
        }
    }

    return network;
}

/* Function to calculate shortest path using Dijkstra's algorithm*/
void calculateShortestPath(const vector<vector<int>>& network, int sourceNode, vector<int>& dist) {
    int numNodes = network.size();
    dist.assign(numNodes, numeric_limits<int>::max());
    dist[sourceNode] = 0;

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push({0, sourceNode});

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        for (int v = 0; v < numNodes; ++v) {
            if (network[u][v] && dist[u] + 1 < dist[v]) {
                dist[v] = dist[u] + 1;
                pq.push({dist[v], v});
            }
        }
    }
}

int calculateTotalSlotsRequired(const std::vector<LightpathRequest>& batch) {
    int totalSlots = 0;
    for (const auto& request : batch) {
        totalSlots += request.getSlotsRequired();
    }
    return totalSlots;
}

std::vector<std::vector<LightpathRequest>> batchFormation(const std::vector<LightpathRequest>& superSet, int batchCapacity) {
    std::vector<std::vector<LightpathRequest>> batches; // Vector to store batches of requests

    std::vector<LightpathRequest> mutableSuperSet(superSet.begin(), superSet.end()); // Create a mutable copy of superSet

    while (!mutableSuperSet.empty()) {
        std::vector<LightpathRequest> currentBatch;
        int remainingCapacity = batchCapacity;

        for (auto it = mutableSuperSet.begin(); it != mutableSuperSet.end();) {
            LightpathRequest& request = *it;

            if (request.getSlotsRequired() <= remainingCapacity) {
                currentBatch.push_back(request);
                remainingCapacity -= request.getSlotsRequired();
                it = mutableSuperSet.erase(it); // Remove the request from mutableSuperSet
            } else {
                ++it;
            }

            if (remainingCapacity <= 0 || it == mutableSuperSet.end()) {
                // Either the batch is full or no more requests can be added to this batch
                batches.push_back(currentBatch);
                break;
            }
        }
    }

    return batches;
}


int main() {
    int numReq, numConnections = 0;
    // Seed the random number generator with current time
    //srand(static_cast<unsigned int>(time(nullptr)));
    //std::srand(static_cast<unsigned int>(std::time(nullptr)));
    std::srand(123); // Replace 123 with the desired fixed seed value

    // define the number of nodes
    const int numNodes = 10;

    // create the netwrk with random connections
    vector<vector<int>> network = createNetwork(numNodes);

    cout << "Number of Nodes = " << numNodes << endl;
    cout << "Total Number of Generated Connections = " << numConnections << endl;
    int totalConnections = (numNodes*(numNodes-1)/2);
    cout << "Total Number of Connections Possible = " << totalConnections << endl;

    int maxSlots = 10; // Adjust this to your specific case
    double mean = maxSlots / 2.0; // Calculate the mean

    // Create a random number generator engine
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a Poisson distribution with the calculated mean
    std::poisson_distribution<int> poisson(mean);


    // Create lightpath requests as objects of the LightpathRequest class
    vector<LightpathRequest> lightpathRequests;
    cout << "Basic Nw : Enter the number of requests to be generated: ";
    cin >> numReq;

    /*Create lightpath requests as object of the LightpathRequest class to save request which are not connected*/
    vector<LightpathRequest> nonConnectedRequests;

    for (int i = 0; i < numReq; ++i) {
        LightpathRequest request;
        request.setId(i + 1);
        request.setSourceNode(rand() % numNodes + 1); // Random source node (1 to numNodes)
        do {
            request.setDestinationNode(rand() % numNodes + 1); // Random destination node (1 to numNodes)
        } while (request.getDestinationNode() == request.getSourceNode()); // Ensure different source and destination

        request.setSlotsRequired(rand() % 10 + 1); // Random number of slots required (1 to 10)
        //request.setSlotsRequired(std::max(1, std::min(poisson(gen), maxSlots))); //using poissons distribution

        // Calculate the number of hops using Dijkstra's algorithm
        vector<int> dist;
        calculateShortestPath(network, request.getSourceNode() - 1, dist);
        //if(dist[request.getDestinationNode() - 1])
        request.setNumHops(dist[request.getDestinationNode() - 1]);
        if (request.getNumHops() > numNodes) {
            request.setConnectionPossible(false);
            nonConnectedRequests.push_back(request);
            // Skip adding this request to lightpathRequests
            continue;
        }
        else {
            request.setConnectionPossible(true);
        }
        lightpathRequests.push_back(request);
    }

    // Print the generated network connections
    cout << "Generated Network Connections:" << endl;
    for (int i = 0; i < numNodes; ++i) {
        for (int j = i + 1; j < numNodes; ++j) {
            if (network[i][j] == 1) {
                cout << "Node " << i + 1 << " is connected to Node " << j + 1 << endl;
            }
        }
    } // */


    // Print the generated lightpath requests
    cout << "\nGenerated Lightpath Requests:" << endl;
    for (const auto& request : lightpathRequests) {
        cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
             << ", Destination: " << request.getDestinationNode() << ", Slots Req:   " << request.getSlotsRequired()
             << ",   No. of Hops: " << request.getNumHops() << " ConctnPosbl?: " << request.isConnectionPossible() << endl;
    } // */

    // Print request which can not be processed
    cout << "\n\nREQUESTS WHICH CAN NOT BE PROCESSED: \n" << endl;
    for (const auto& request : nonConnectedRequests) {
        cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
             << ", Destination: " << request.getDestinationNode() << ", Slots Req:   " << request.getSlotsRequired()
             << ",   No. of Hops: " << request.getNumHops() << " ConctnPosbl?: " << request.isConnectionPossible() << endl;
    }

    // Calculate maximum and minimum numHops
    //int hmax = numeric_limits<int>::min();
    //int hmin = numeric_limits<int>::max();


    int hmax = 1;
    int hmin = numNodes;

    for (const auto& request : lightpathRequests) {
        int numHops = request.getNumHops();
        if (numHops > hmax) {
            hmax = numHops;
        }
        if (numHops < hmin) {
            hmin = numHops;
        }
    }

    cout << "hmin = " << hmin << endl;
    cout << "hmax = " << hmax << endl;

    // Calculate the number of variables to create
    int numSets = hmax - hmin + 1;

    // Create an array of vectors to store the sets
    vector<LightpathRequest> sets[numSets]; /* THESE ARE THE CLASSES WHICH WILL HOLD REQUESTS WITH  Particular number of Hops in each class/set*/

    // Copy members of lightpathRequests to respective sets
    for (const auto& request : lightpathRequests) {
        int numHops = request.getNumHops();
        sets[numHops - hmin].push_back(request);
    }


    /*// Print the sets
    for (int i = 0; i < numSets; ++i) {
        cout << "Set " << i + 1 << ":" << endl;
        for (const auto& request : sets[i]) {
            cout << "ID: " << request.getId() << ", Source : " << request.getSourceNode()
                 << ", Destination : " << request.getDestinationNode() << ", Slots Required: " << request.getSlotsRequired()
                 << ", Number of Hops: " << request.getNumHops() << endl;
        }
        cout << endl;
    } // */

    // Print  sets SIZE
    for (int i = 0; i < numSets; ++i) {
        cout << "Set " << i+1 << " Size =" << sets[i].size() << endl;;
    }

    // After populating the sets, sort the requests within each set based on slotsRequired
    for (int i = 0; i < numSets; ++i) {
        std::sort(sets[i].begin(), sets[i].end(), [](const LightpathRequest& a, const LightpathRequest& b) {
            return a.getSlotsRequired() < b.getSlotsRequired();
        });
    }

    /*// Print the sorted sets
    for (int i = 0; i < numSets; ++i) {
        cout << "Set " << i + 1 << " (sorted by slotsRequired):" << endl;
        for (const auto& request : sets[i]) {
            cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
                    << ", Destination: " << request.getDestinationNode() << ", Slots Required: " << request.getSlotsRequired()
                    << ", Number of Hops: " << request.getNumHops() << endl;
        }
        cout << endl;
    } // */

    vector<LightpathRequest> superSet;

    // Copy elements from sets[numSets] to sets[0] into superSet
    for (int i = numSets - 1; i >= 0; --i) {
        superSet.insert(superSet.end(), sets[i].begin(), sets[i].end());
    }

    int supSetCount =0;
    //Printing Super Set
    cout << "\n\n\n Superset : \n" << endl;
    for (const auto& request : superSet) {
        cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
             << ", Destination: " << request.getDestinationNode() << ", Slots Req:   " << request.getSlotsRequired()
             << ",   No. of Hops: " << request.getNumHops() << " ConctnPosbl?: " << request.isConnectionPossible() << endl;
            supSetCount++;
    } // */
    /*for (const auto& request : superSet) {
        supSetCount++;f
    }//*/
    cout << "Number of elements in superSet (total no of requests) = " << supSetCount << endl;

    int batchCapacity = 40;
    // Call batchFormation function and save the batches in the "batches" vector
    std::vector<std::vector<LightpathRequest>> batches = batchFormation(superSet, batchCapacity);

    // Print the batches
    for (int i = 0; i < batches.size(); ++i) {
        std::cout << "Batch " << i + 1 << " (Total Slots Required: " << calculateTotalSlotsRequired(batches[i]) << "):" << std::endl;
        for (const auto& request : batches[i]) {
            std::cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
                 << ", Destination: " << request.getDestinationNode() << ", Slots Req: " << request.getSlotsRequired()
                 << ", No. of Hops: " << request.getNumHops() << " ConctnPosbl?: " << request.isConnectionPossible() << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
