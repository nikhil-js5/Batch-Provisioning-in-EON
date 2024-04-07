#include <iostream>
#include <vector>
#include <array>
#include <queue>
#include <limits>
#include <functional>
#include <algorithm>
#include <ctime>   // for srand()
using namespace std;

const int num_nodes = 14;
int MAX_SLOTS = 320; // Initialize network slots as all free
int slot_size = 30;
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
public:
    // Public member functions and data members
    int id;
    int sourceNode;
    int destinationNode;
    int slotsRequired;
    int numHops;
    bool connectionPossible;

    LightpathRequest() : id(0), sourceNode(0), destinationNode(0),
                         slotsRequired(0), numHops(0), connectionPossible(true) {}

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
    for (int i = 0; i < numRequests; ++i) {
        LightpathRequest req;
        req.id = i + 1;
        req.sourceNode = rand() % num_nodes + 1; // Generate random Src node
        do {
            req.destinationNode = rand() % num_nodes + 1;  // Generate random destination node
        } while (req.destinationNode == req.sourceNode);  // Ensure source and destination are not the same


        req.slotsRequired = rand() % slot_size + 1; // Random slots required, max 320
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
// SPFF Algorithm - Assigns the first available contiguous block of spectrum slots to the request
bool SPFF_Algorithm(int sourceNode, int destinationNode, int slotsRequired, AdjacencyMatrix& network, AdjacencyMatrix& slots) {
    // Find the shortest path between source and destination
    std::vector<int> path = findShortestPath(sourceNode, destinationNode, network);
    /*for(int i=0;i < path.size() ; i++){
        std::cout << path[i] +1 << " ->";
    }
    cout << endl; // */

    // Check for available contiguous slots on the path
    if (!path.empty() && checkContiguousSlots(path, slotsRequired, network, slots)) {
        // Assign the slots if contiguous slots are available
        return assignSlotsToRequest(path, slotsRequired, network, slots);
    }

    // If there's no path or not enough contiguous slots, return false indicating the request is blocked
    return false;
}


int main() {
    // Seeding the random number generator with current time
    //std::srand(static_cast<unsigned int>(std::time(nullptr)));
   std::srand(123); // fixed seed value

    int numRequests = 100;
    for(numRequests = 100; numRequests<2001; numRequests+=100){
            cout << numRequests << " Request \n";
    AdjacencyMatrix network = {{
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
    }};
    /*AdjacencyMatrix network = {
        {0,1,1,1,0},
        {1,0,1,0,0},
        {1,1,0,1,1},
        {1,0,1,0,0},
        {0,0,1,0,0},
    };*/

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

    int hmax = 1;
    int hmin = num_nodes;

    for (const auto& request : lightpathRequests) {
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
        superSet.insert(superSet.end(), sets[i].begin(), sets[i].end());
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
    std::vector<std::vector<LightpathRequest>> batches = batchFormation(superSet, MAX_SLOTS);

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


    //----------
    //std::vector<std::vector<int>> slots(num_nodes, std::vector<int>(MAX_SLOTS, 0));
    //AdjacencyMatrix slots(num_nodes, std::vector<int>(MAX_SLOTS, 0));
    AdjacencyMatrix slots = initializeSlotMatrix(network, MAX_SLOTS);

    int totalBlocked = 0, ServedRequests =0;
    int totalRequests = lightpathRequests.size();

    // Check each lightpath request for available slots and blocking
    for (int i = 0; i < batches.size(); ++i){
        for (const auto& request : batches[i]) {
            bool slotsAssigned = SPFF_Algorithm(request.getSourceNode() -1, request.getDestinationNode() -1, request.getSlotsRequired(), network, slots);
            //printSlots(network, slots);
            if (slotsAssigned) {
                // Slots are available, request is not blocked
                ++ServedRequests;
               // cout << "batch " << i << " served request : " << request.id << endl;
            } else {
                // No slots available, increment blocked requests count
                ++totalBlocked;
                //cout << "batch " << i <<  " rejected request : " << request.id  << endl ;
            }
        }
    }
    // Compute and display the blocking probability
    double blockingProbability = static_cast<double>(totalBlocked) / totalRequests;
    //std::cout << "Served requests = " << ServedRequests << " : Blocked Requests = " << totalBlocked << std::endl;
    std::cout << "\nBlocking Probability: WITH Batch " << blockingProbability << std::endl;

    AdjacencyMatrix slots_1 = initializeSlotMatrix(network, MAX_SLOTS);
    totalBlocked = 0; ServedRequests =0;
    totalRequests = lightpathRequests.size();
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
    double blockingProbability1 = static_cast<double>(totalBlocked) / totalRequests;
    //std::cout << "Served requests " << ServedRequests << " : Blocked Requests = " << totalBlocked << std::endl;
    std::cout << "Blocking Probability without batch " << blockingProbability1 << std::endl;
    //std::cout << "Blocking Probability: WITH Batch " << blockingProbability << std::endl;



    }
    return 0;
}
