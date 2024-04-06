#include <iostream>
#include <vector>
#include <array>
#include <queue>
#include <limits>
#include <functional>
#include <algorithm>

const int num_nodes = 14;
using AdjacencyMatrix = std::vector<std::vector<int>>;
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

    bool operator<(const LightpathRequest& rhs) const {
        // Prioritize by least number of hops and then by fewest slots required
        return (numHops < rhs.numHops) || (numHops == rhs.numHops && slotsRequired < rhs.slotsRequired);
    }
};

std::vector<LightpathRequest> generateLightpathRequests(int num_nodes, int numRequests);
std::vector<std::vector<LightpathRequest>> batchFormation(const std::vector<LightpathRequest>& superSet, int batchCapacity);

// Function to check for available contiguous slots on the path
bool checkContiguousSlots(const std::vector<int>& path, int slotsRequired, const std::vector<std::vector<int>>& network) {
    // Iterate over each link in the path
    for (const auto& link : path) {
        int freeSlots = 0;
        for (int slot : network[link]) {
            freeSlots = (slot == 0) ? (freeSlots + 1) : 0;
            if (freeSlots == slotsRequired) {
                return true;
            }
        }
    }
    return false;
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



std::vector<LightpathRequest> generateLightpathRequests(int num_nodes, int numRequests) {
    std::vector<LightpathRequest> requests;
    for (int i = 0; i < numRequests; ++i) {
        LightpathRequest req;
        req.id = i + 1;
        req.sourceNode = rand() % num_nodes + 1;
        req.destinationNode = rand() % num_nodes + 1;
        req.slotsRequired = rand() % 50 + 1; // Random slots required, max 320
        req.numHops = rand() % 10 + 1;  // Random number of hops, for example purposes
        requests.push_back(req);
    }
    return requests;
}

///
bool assignSlotsToRequest(const std::vector<int>& path, int slotsRequired, std::vector<std::vector<int>>& network) {
    // For simplicity, we're assuming the network is represented as a 2D vector of int,
    // where each element represents a slot that can be 0 (free) or 1 (occupied).
    for (const auto& link : path) {
        int freeSlots = 0;
        int startIndex = -1;
        // Find the first fit of contiguous slots for the request
        for (int i = 0; i < network[link].size(); ++i) {
            if (network[link][i] == 0) {  // Slot is free
                if (startIndex == -1) {
                    startIndex = i;  // Start of a possible block of free slots
                }
                freeSlots++;
                if (freeSlots == slotsRequired) {
                    // We found a block of free slots, mark them as occupied
                    for (int j = startIndex; j < startIndex + slotsRequired; ++j) {
                        network[link][j] = 1;  // Mark as occupied
                    }
                    return true;  // The slots were successfully assigned
                }
            } else {
                freeSlots = 0;  // Reset the counter if we hit an occupied slot
                startIndex = -1;  // Reset the start index
            }
        }
    }
    return false;  // Not enough contiguous free slots were found for this request
}


// SPFF Algorithm - Assigns the first available contiguous block of spectrum slots to the request
bool SPFF_Algorithm(int sourceNode, int destinationNode, int slotsRequired, AdjacencyMatrix& network) {
    // Find the shortest path between source and destination
    std::vector<int> path = findShortestPath(sourceNode, destinationNode, network);

    // Check for available contiguous slots on the path
    if (!path.empty()) {
        if (checkContiguousSlots(path, slotsRequired, network)) {
            // Assign the slots
            // Note: Here we need to actually assign the slots to the request in the 'network' data structure.
            // This implementation assumes that a function will assign the slots.
            assignSlotsToRequest(path, slotsRequired, network);

            return true;
        }
    }
    return false; // Request is blocked if there are no contiguous slots
}

int main() {
    int numRequests = 10000;
   // AdjacencyMatrix network = {};
    const AdjacencyMatrix network = {{
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

    // Generate lightpath requests
    std::vector<LightpathRequest> lightpathRequests = generateLightpathRequests(num_nodes, numRequests);

    // Initialize network slots as all free
    int MAX_SLOTS = 320;
    std::vector<std::vector<int>> slots(num_nodes, std::vector<int>(MAX_SLOTS, 0));


    int totalBlocked = 0;
    int totalRequests = lightpathRequests.size();

    // Check each lightpath request for available slots and blocking
    for (const auto& request : lightpathRequests) {
        std::vector<int> path = findShortestPath(request.getSourceNode() - 1, request.getDestinationNode() - 1, network);

        // If there's a valid path, check for contiguous slots
        if (!path.empty() && checkContiguousSlots(path, request.getSlotsRequired(), network)) {
            // Slots are available, request is not blocked
        } else {
            // No slots available, increment blocked requests count
            ++totalBlocked;
        }
    }

    // Compute and display the blocking probability
    double blockingProbability = static_cast<double>(totalBlocked) / totalRequests;
    std::cout << "Blocking Probability: " << blockingProbability << std::endl;

    return 0;
}
