#include <iostream>
#include <vector>
#include <array>
#include <queue>
#include <limits>
#include <functional>
#include <algorithm>
#include <ctime>   // for srand()
#include <fstream> // for saving data into external file
#include <random>

using namespace std;

const int num_nodes = 28;
int SLOTS_AVAILABLE = 320;

int slot_size = 20;
using AdjacencyMatrix = std::vector<std::vector<int>>;
class LightpathRequest;


AdjacencyMatrix initializeSlotMatrix(const AdjacencyMatrix& network, int max_slots) {
    AdjacencyMatrix slots(network.size(), std::vector<int>(network.size(), -1));
    for (int i = 0; i < network.size(); i++) {
        for (int j = 0; j < network[i].size(); j++) {
            if (network[i][j] != 0) {
                slots[i][j] = max_slots;  // Set the maximum number of slots for each link
            }
        }
    }
    return slots;
}


std::vector<LightpathRequest> generateLightpathRequests(int num_nodes, int numRequests);
std::vector<std::vector<LightpathRequest>> batchFormation(const std::vector<LightpathRequest>& superSet, int batchCapacity);

// Function to check for available contiguous slots on the path
bool checkContiguousSlots(const std::vector<int>& path, int slotsRequired, const std::vector<std::vector<int>>& network, AdjacencyMatrix& slots);
//Dijkstra Shortest path algorithm
std::vector<int> findShortestPath(int source, int target, const AdjacencyMatrix& network);
// Assigning the slots
bool assignSlotsToRequest(const std::vector<int>& path, int slotsRequired, std::vector<std::vector<int>>& network, AdjacencyMatrix& slots) ;
// SPFF Algorithm - Assigns the first available contiguous block of spectrum slots to the request
bool SPFF_Algorithm(int sourceNode, int destinationNode, int slotsRequired, AdjacencyMatrix& network, AdjacencyMatrix& slots) ;

//
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

           double getReqTime() const {
        return req_time;
    }

    void setReqTime(double newReqTime) {
        req_time = newReqTime;
    }

public:
    // Public member functions and data members
    int id;
    int sourceNode;
    int destinationNode;
    int slotsRequired;
    int numHops;
    bool connectionPossible;
    double req_time;  // Arrival time of the request

    LightpathRequest() : id(0), sourceNode(0), destinationNode(0),
                         slotsRequired(0), numHops(0), connectionPossible(true), req_time(0.0) {}


    /*bool operator<(const LightpathRequest& rhs) const {
        // Prioritize by least number of hops and then by fewest slots required
        return (numHops < rhs.numHops) || (numHops == rhs.numHops && slotsRequired < rhs.slotsRequired);
    }*/
};


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


void printSlots(const AdjacencyMatrix& network, const AdjacencyMatrix& slots) {
    cout << "\nSlots status:" << endl;
    for (int i = 0; i < network.size(); ++i) {
        for (int j = 0; j < network[i].size(); ++j) {
            if (network[i][j] != 0) {  // There is a link between node i and node j
                cout << "Edge " << i + 1 << " -> " << j + 1 << ": ";
                if (slots[i][j] >= 0) {  // Check if slots are initialized for this edge
                    cout << slots[i][j] << " slots available" << endl;
                } else {
                    cout << "No slots available (unconnected)" << endl;
                }
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

std::vector<LightpathRequest> generateLightpathRequests(int num_nodes, int numRequests) {
    std::vector<LightpathRequest> requests;

    std::random_device rd;
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::exponential_distribution<> d(1.0); // Mean arrival rate (lambda = 1)

    double currentTime = 0.0;


    for (int i = 0; i < numRequests; ++i) {
        LightpathRequest req;
        req.id = i + 1;
        req.sourceNode = rand() % num_nodes + 1; // Generate random Src node
        do {
            req.destinationNode = rand() % num_nodes + 1;  // Generate random destination node
        } while (req.destinationNode == req.sourceNode);  // Ensure source and destination are not the same


        req.slotsRequired = rand() % slot_size + 1; // Random slots required, max 320

        double interArrivalTime = d(gen); // Generate the inter-arrival time
        currentTime += interArrivalTime; // Update the current time by adding the inter-arrival time
        req.setReqTime(currentTime); // Set the arrival time of the request


       // req.numHops = rand() % 10 + 1;  // Random number of hops, for example purposes
        requests.push_back(req);
    }
    // Print the generated lightpath requests
    /*cout << "\nGenerated Lightpath Requests:" << endl;
    for (const auto& request : requests) {
        cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
             << ", Destination: " << request.getDestinationNode() << ", Slots Req:   " << request.getSlotsRequired()
             << ",   No. of Hops: " << request.getNumHops()  << endl;
    } // */
    return requests;
}


std::vector<int> findShortestPath(int source, int target, const AdjacencyMatrix& network) {
    std::vector<int> dist(num_nodes, std::numeric_limits<int>::max());
    std::vector<int> prev(num_nodes, -1);
    std::vector<bool> visited(num_nodes, false);

    using pii = std::pair<int, int>;
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> queue;

    dist[source] = 0;
    queue.push({0, source});

    while (!queue.empty()) {
        int distance = queue.top().first;
        int current = queue.top().second;
        queue.pop();

        if (visited[current]) continue;
        visited[current] = true;

        for (int i = 0; i < num_nodes; ++i) {
            if (network[current][i] && !visited[i]) {
                int newDist = distance + 1;  // Assuming each edge has a weight of 1
                if (newDist < dist[i]) {
                    dist[i] = newDist;
                    prev[i] = current;
                    queue.push({newDist, i});
                }
            }
        }
    }

    // Reconstruct the shortest path from source to target
    std::vector<int> path;
    if (dist[target] == std::numeric_limits<int>::max()) {
        return path; // No path exists
    }
    for (int at = target; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    std::reverse(path.begin(), path.end());
    return path;
}


bool checkContiguousSlots(const std::vector<int>& path, int slotsRequired, const AdjacencyMatrix& network, AdjacencyMatrix& slots) {
    // Assume all edges in the path have at least the required slots available
    for (size_t i = 0; i < path.size() - 1; i++) {
        int startNode = path[i];
        int endNode = path[i + 1];
        if (slots[startNode][endNode] < slotsRequired) {
            return false;  // Not enough slots available on this edge
        }
    }
    return true;  // Sufficient contiguous slots available on all edges
}

bool assignSlotsToRequest(const std::vector<int>& path, int slotsRequired, AdjacencyMatrix& network, AdjacencyMatrix& slots) {
    if (checkContiguousSlots(path, slotsRequired, network, slots)) {
        for (size_t i = 0; i < path.size() - 1; i++) {
            int startNode = path[i];
            int endNode = path[i + 1];
            slots[startNode][endNode] -= slotsRequired;  // Reduce available slots by the amount required
        }
        return true;
    }
    return false;  // Not enough contiguous slots were found for this request
}





// SPFF Algorithm - Assigns the first available contiguous block of spectrum slots to the request
/*bool SPFF_Algorithm(int sourceNode, int destinationNode, int slotsRequired, AdjacencyMatrix& network, AdjacencyMatrix& slots) {
    // Find the shortest path between source and destination
    std::vector<int> path = findShortestPath(sourceNode, destinationNode, network);

    // Check for available contiguous slots on the path
    if (!path.empty()) {
        if (checkContiguousSlots(path, slotsRequired, network, slots)) {
            // Assign the slots
            // Note: Here we need to actually assign the slots to the request in the 'network' data structure.
            // This implementation assumes that a function will assign the slots.
            assignSlotsToRequest(path, slotsRequired, network, slots);

            return true;
        }
    }
    return false; // Request is blocked if there are no contiguous slots
}*/
/*
// SPFF Algorithm - Assigns the first available contiguous block of spectrum slots to the request
bool SPFF_Algorithm(int sourceNode, int destinationNode, int slotsRequired, AdjacencyMatrix& network, AdjacencyMatrix& slots) {
    // Find the shortest path between source and destination
    std::vector<int> path = findShortestPath(sourceNode, destinationNode, network);
    /*for(int i=0;i < path.size() ; i++){
        std::cout << path[i] +1 << " ->";
    }
    cout << endl; //

    // Check for available contiguous slots on the path
    if (!path.empty() && checkContiguousSlots(path, slotsRequired, network, slots)) {
        // Assign the slots if contiguous slots are available
        return assignSlotsToRequest(path, slotsRequired, network, slots);
    }

    // If there's no path or not enough contiguous slots, return false indicating the request is blocked
    return false;
}*/

bool SPFF_Algorithm(int sourceNode, int destinationNode, int slotsRequired, AdjacencyMatrix& network, AdjacencyMatrix& slots) {
    // Find the shortest path between source and destination
    std::vector<int> path = findShortestPath(sourceNode, destinationNode, network);

    //if(request.)

    if (!path.empty() && checkContiguousSlots(path, slotsRequired, network, slots)) {
        return assignSlotsToRequest(path, slotsRequired, network, slots);
    }

    // If the first path is not usable, attempt to find an alternative path
    // Increase the weights of the edges used in the first path to discourage their selection
    if (!path.empty()) {
        for (size_t i = 0; i < path.size() - 1; i++) {
            int startNode = path[i];
            int endNode = path[i + 1];
            network[startNode][endNode] += 10000; // Adding a large value to the weight to make it non-preferable
            network[endNode][startNode] += 10000; // Ensure bidirectional paths are equally adjusted if undirected
        }

        // Try finding another shortest path with updated network weights
        std::vector<int> secondPath = findShortestPath(sourceNode, destinationNode, network);
        if (!secondPath.empty() && checkContiguousSlots(secondPath, slotsRequired, network, slots)) {
            bool result = assignSlotsToRequest(secondPath, slotsRequired, network, slots);
            // Restore original weights
            for (size_t i = 0; i < path.size() - 1; i++) {
                int startNode = path[i];
                int endNode = path[i + 1];
                network[startNode][endNode] -= 10000;
                network[endNode][startNode] -= 10000;
            }
            return result;
        }

        // Restore original weights regardless of second path success
        for (size_t i = 0; i < path.size() - 1; i++) {
            int startNode = path[i];
            int endNode = path[i + 1];
            network[startNode][endNode] -= 9999;
            network[endNode][startNode] -= 9999;
        }
    }//*/

    // If no viable path was found or not enough contiguous slots, return false
    return false;
}

int main() {
    // Seeding the random number generator with current time
    //std::srand(static_cast<unsigned int>(std::time(nullptr)));
   std::srand(50); // fixed seed value



    std::ofstream outFile("blocking_probabilities.txt");
    if (!outFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return 1;
    }

       /* AdjacencyMatrix network = {{ //: NSF network
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1},
        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}
    }};  // */

      /*  AdjacencyMatrix network = {{ //Indian Network
        {0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0},
    }};*/

        AdjacencyMatrix network = {{ //European Network
        {0, 0, 1, 0, 0, 1, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 1, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0, 0, 1,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 1, 0, 0, 1,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 1, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 1, 1, 0,       0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0,       1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0, 0, 1,       0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 0, 0, 0, 1, 0,       0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0,       0, 1, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},

        {0, 0, 0, 0, 0, 0, 1, 0, 0, 0,       1, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1,       0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0,       0, 1, 1, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       1, 0, 0, 0, 0, 1, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 1, 1, 0, 1, 0, 0, 0,  1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 1, 0, 0, 1, 0, 1, 0, 0,  0, 1, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 1, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 1, 1, 0, 0, 0, 0, 0},

        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  1, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  1, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 1, 0, 0, 0, 1, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 1, 0, 1, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0, 0, 0},
    }};
    /*AdjacencyMatrix network = {
        {0,1,1,1,0},
        {1,0,1,0,0},
        {1,1,0,1,1},
        {1,0,1,0,0},
        {0,0,1,0,0},
    };*/

    int numRequests;
    for(numRequests = 4000; numRequests<=4000; numRequests+=100){
            cout << numRequests << " Requests \n";


    // Generate lightpath requests
    std::vector<LightpathRequest> lightpathRequests = generateLightpathRequests(num_nodes, numRequests);

    for (auto& request : lightpathRequests) {
        std::vector<int> path = findShortestPath(request.getSourceNode() - 1, request.getDestinationNode() - 1, network);
        if (path.empty()) {
            // If no path is found, set numHops to a default or error value, e.g., -1
            request.setNumHops(-1);
        } else {
            // Set the number of hops; subtract one because the path includes both the start and end nodes
            request.setNumHops(path.size() - 1);
        }
    }


   /* // Print the generated lightpath requests
    cout << "\nGenerated Lightpath Requests:" << endl;
    for (const auto& request : lightpathRequests) {
        cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
             << ", Destination: " << request.getDestinationNode() << ", Slots Req:   " << request.getSlotsRequired()
             << ",   No. of Hops: " << request.getNumHops()  << endl;
    } // */


    //////// WITHOUT BATCH PROCESSING
    AdjacencyMatrix slots_1 = initializeSlotMatrix(network, SLOTS_AVAILABLE);
    int totalBlocked = 0, ServedRequests =0;
    int totalRequests = lightpathRequests.size();
    // Check each lightpath request for available slots and blocking
    for (const auto& request : lightpathRequests) {
        bool slotsAssigned = SPFF_Algorithm(request.getSourceNode() -1, request.getDestinationNode() -1, request.getSlotsRequired(), network, slots_1);
        //printSlots(network, slots);
        if (slotsAssigned) {
            // Slots are available, request is not blocked
            ++ServedRequests;
          //  cout << "served request WITH BATCH: " << request.id << endl;
        } else {
            // No slots available, increment blocked requests count
            ++totalBlocked;
          //  cout << " rejected request : WITH BATCH " << request.id  << endl ;
        }
    }
    // Compute and display the blocking probability
    double blockingProbability = static_cast<double>(totalBlocked) / totalRequests;
    //std::cout << "Served requests " << ServedRequests << " : Blocked Requests = " << totalBlocked << std::endl;
    std::cout << "Blocking Probability without batch: " << blockingProbability ;
    //std::cout << "Blocking Probability: WITH Batch " << blockingProbability << std::endl;
    ///outFile << blockingProbability << std::endl ; // without Batch



    int MAX_SLOTS ; // Initialize network slots as all free
    //std::vector<int> MAX_REQUESTS = {1, 10, 20, 50, 100, 200, 300, 400, 500, 700,900, 1000, 1500,2000};
    //std::vector<int> MAX_REQUESTS = {1, 10, 20, 50, 100, 200, 300, 400, 500, 700,900, 1000, 1500,2000};


    //for(int MAX_SLOTS : MAX_REQUESTS){
    ///for(MAX_SLOTS = 200; MAX_SLOTS <= numRequests; MAX_SLOTS +=20){
    for(MAX_SLOTS = 100; MAX_SLOTS <= 4000; MAX_SLOTS +=100){
                ////// BATCH PROCESSING OF REQUESTS
            int totalBlocked_B = 0, ServedRequests_B =0;
            //totalRequests = lightpathRequests.size();
            double blockingProbability_B;

    //std::vector<std::vector<int>>    slots(num_nodes, std::vector<int>(MAX_SLOTS, 0));
    //AdjacencyMatrix slots(num_nodes, std::vector<int>(MAX_SLOTS, 0));
    AdjacencyMatrix slots = initializeSlotMatrix(network, SLOTS_AVAILABLE);

    for(int i=0;i<totalRequests;i+=MAX_SLOTS){
    // Define the end of the current batch, making sure not to exceed the vector's size
    //size_t end = std::min(lightpathRequests.size(), i + MAX_SLOTS);
    size_t end = std::min(lightpathRequests.size(), static_cast<size_t>(i + MAX_SLOTS));


    // Create a batch using iterators from the current segment of the vector
    std::vector<LightpathRequest> currentBatch(lightpathRequests.begin() + i, lightpathRequests.begin() + end);



    int hmax = 1;
    int hmin = num_nodes;

    for (const auto& request : currentBatch) {
        int numHops = request.getNumHops();
        if (numHops > hmax) {
            hmax = numHops;
        }
        if (numHops < hmin) {
            hmin = numHops;
        }
    }

    //cout << "hmin = " << hmin << endl;
    //cout << "hmax = " << hmax << endl;

    // Calculate the number of variables to create
    int numSets = hmax - hmin + 1;

     // Create an array of vectors to store the sets
    vector<LightpathRequest> sets[numSets]; /* THESE ARE THE CLASSES WHICH WILL HOLD REQUESTS WITH  Particular number of Hops in each class/set*/

    // Copy members of lightpathRequests to respective sets
    for (const auto& request : currentBatch) {
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

   /* // Print  sets SIZE
    for (int i = 0; i < numSets; ++i) {
        cout << "Set " << i+1 << " Size =" << sets[i].size() << endl;;
    } // */

    //----------
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
        superSet.insert(superSet.begin(), sets[i].begin(), sets[i].end());
    }

    int supSetCount =0;
    //Printing Super Set
  //  cout << "\n\n\n Superset : \n" << endl;
    for (const auto& request : superSet) {
       /* cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
             << ", Destination: " << request.getDestinationNode() << ", Slots Req:   " << request.getSlotsRequired()
             << ",   No. of Hops: " << request.getNumHops() << " ConctnPosbl?: " << request.isConnectionPossible() << endl;*/
            supSetCount++;
    } // */
    /*for (const auto& request : superSet) {
        supSetCount++;f
    }
    cout << "Number of elements in superSet (total no of requests) = " << supSetCount << endl;//*/

    // Call batchFormation function and save the batches in the "batches" vector
    std::vector<std::vector<LightpathRequest>> batches = batchFormation(superSet, SLOTS_AVAILABLE);

     // Print the batches
   /* for (int i = 0; i < batches.size(); ++i) {
        std::cout << "Batch " << i + 1 << " (Total Slots Required: " << calculateTotalSlotsRequired(batches[i]) << "):" << std::endl;
        for (const auto& request : batches[i]) {
            std::cout << "ID: " << request.getId() << ", Source: " << request.getSourceNode()
                 << ", Destination: " << request.getDestinationNode() << ", Slots Req: " << request.getSlotsRequired()
                 << ", No. of Hops: " << request.getNumHops() << " ConctnPosbl?: " << request.isConnectionPossible() << std::endl;
        }
        std::cout << std::endl;
    } // */





    // Check each lightpath request for available slots and blocking
    for (int i = 0; i < batches.size(); ++i){
       // for (int j = 0; i < batches.size(); ++i)
        for (const auto& request : batches[i]) {
            bool slotsAssigned = SPFF_Algorithm(request.getSourceNode() -1, request.getDestinationNode() -1, request.getSlotsRequired(), network, slots);
            //printSlots(network, slots);
            if (slotsAssigned) {
                // Slots are available, request is not blocked
                ++ServedRequests_B;
               // cout << "batch " << i << " served request : " << request.id << endl;
            } else {
                // No slots available, increment blocked requests count
                ++totalBlocked_B;
                //cout << "batch " << i <<  " rejected request : " << request.id  << endl ;
            }
        }
    }

    }
    // Compute and display the blocking probability
    blockingProbability_B = static_cast<double>(totalBlocked_B) / totalRequests;
    //std::cout << "Served requests = " << ServedRequests << " : Blocked Requests = " << totalBlocked << std::endl;
    std::cout << "\n For Batch size ::: " << MAX_SLOTS << " ::: BP: WITH Batch:: " << blockingProbability_B << std::endl ;

    // Write the computed probabilities to the file
    //outFile << "without   : " << blockingProbability1 << std::endl;
    //outFile << blockingProbability << std::endl ; // without Batch
    outFile << blockingProbability_B << std::endl ; // with Batch


    }
    }
    //}
    outFile.close();
    return 0;
}
