#include <iostream>
#include <vector>
#include <map>

class OpticalNode {
public:
    OpticalNode(int id) : id(id) {}

    int getId() const {
        return id;
    }

private:
    int id;
};

class OpticalConnection {
public:
    OpticalConnection(int source, int destination, int bandwidth) :
        source(source), destination(destination), bandwidth(bandwidth) {}

    int getSource() const {
        return source;
    }

    int getDestination() const {
        return destination;
    }

    int getBandwidth() const {
        return bandwidth;
    }

private:
    int source;
    int destination;
    int bandwidth;
};

class OpticalLink {
public:
    // Default constructor with reasonable default values
    OpticalLink() : source(-1), destination(-1), numChannels(0), allocatedChannels(0) {}

    OpticalLink(int source, int destination, int numChannels) :
        source(source), destination(destination), numChannels(numChannels), allocatedChannels(0) {
        // Initialize spectrum allocation
        spectrumAllocation.resize(numChannels, false); // false to indicate unused channel
    }

    int getSource() const {
        return source;
    }

    int getDestination() const {
        return destination;
    }

    int getNumChannels() const {
        return numChannels;
    }

    int getAllocatedChannels() const {
        return allocatedChannels;
    }

    bool allocateSpectrum(int bandwidth) {
        // Simple FIRST FIT spectrum allocation
        int requiredChannels = (bandwidth + 9) / 10; // Each channel has 10 units of bandwidth
        for (int i = 0; i <= numChannels - requiredChannels; ++i) {
            bool canAllocate = true;
            for (int j = i; j < i + requiredChannels; ++j) {
                if (spectrumAllocation[j]) {
                    canAllocate = false;
                    break;
                }
            }
            if (canAllocate) {
                for (int j = i; j < i + requiredChannels; ++j) {
                    spectrumAllocation[j] = true;
                }
                allocatedChannels += requiredChannels;
                return true;
            }
        }
        return false;
    }

private:
    int source;
    int destination;
    int numChannels;
    int allocatedChannels;
    std::vector<bool> spectrumAllocation;
};

int main() {
    // Create optical network "nodes"
    std::vector<OpticalNode> nodes;
    for (int i = 1; i <= 5; ++i) {
        nodes.push_back(OpticalNode(i));
    }

    // Create optical connections
    std::vector<OpticalConnection> connections;
    connections.push_back(OpticalConnection(1, 2, 10));
    connections.push_back(OpticalConnection(1, 2, 25));
    connections.push_back(OpticalConnection(1, 2, 20));
    connections.push_back(OpticalConnection(2, 3, 7));
    connections.push_back(OpticalConnection(3, 4, 120));
    connections.push_back(OpticalConnection(4, 5, 8));
    connections.push_back(OpticalConnection(1, 4, 15));
    connections.push_back(OpticalConnection(2, 5, 9));
    connections.push_back(OpticalConnection(3, 5, 11));

    // Create optical links between nodes with spectrum allocation
    std::map<std::pair<int, int>, OpticalLink> opticalLinks; /* Source and
    Destination are Key and OpticalLink is the value for the MAP*/
    opticalLinks[{1, 2}] = OpticalLink(1, 2, 5);  // 100 channels
    opticalLinks[{2, 3}] = OpticalLink(2, 3, 80);   // 80 channels
    opticalLinks[{3, 4}] = OpticalLink(3, 4, 12);  // 120 channels
    opticalLinks[{4, 5}] = OpticalLink(4, 5, 60);   // 60 channels
    opticalLinks[{1, 4}] = OpticalLink(1, 4, 150);  // 150 channels
    opticalLinks[{2, 5}] = OpticalLink(2, 5, 90);   // 90 channels
    opticalLinks[{3, 5}] = OpticalLink(3, 5, 110);  // 110 channels

    // Allocate spectrum for connections
    for (const OpticalConnection& connection : connections) {
        int source = connection.getSource();
        int destination = connection.getDestination();
        int bandwidth = connection.getBandwidth();

        if (opticalLinks.count({source, destination}) > 0) {
        /*'.count({source, destination}): This is a method call on the opticalLinks map.
        It checks how many times a specific key, represented as the pair {source, destination},
        appears in the map. The count method returns the number of occurrences, which
        will be either 0 (if the key is not found) or 1 (if the key is found, as map keys are
        unique). the line of code is saying: "If the key {source, destination} exists at least
        once in the opticalLinks map, then execute the code within the following block."*/
            OpticalLink& link = opticalLinks[{source, destination}];
            if (link.allocateSpectrum(bandwidth)) {
                std::cout << "Connection from Node " << source << " to Node " << destination
                          << " with bandwidth " << bandwidth << " allocated successfully. Allocated Channels of link = 0" << link.getAllocatedChannels() <<  std::endl;
            } else {
                std::cout << "Connection from Node " << source << " to Node " << destination
                          << " with bandwidth " << bandwidth << " cannot be allocated due to insufficient spectrum. Available Channels in link are : " << (link.getNumChannels() - link.getAllocatedChannels()) << std::endl;
            }
        } else {
            std::cout << "No optical link between Node " << source << " and Node " << destination << "." << std::endl;
        }
    }

    return 0;
}
