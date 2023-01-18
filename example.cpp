#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void generateRoutes(vector<int> &path, int start) {
    if (path.size() == 5) {
        for (int i = 0; i < path.size(); i++) {
            cout << path[i] << " ";
        }
        cout << endl;
        return;
    }

    for (int i = 0; i < 5; i++) {
        if (find(path.begin(), path.end(), i) == path.end()) {
            path.push_back(i);
            generateRoutes(path, i);
            path.pop_back();
        }
    }
}

int main() {
    vector<int> path;
    path.push_back(0);
    generateRoutes(path, 0);
    return 0;
}