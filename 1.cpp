#include <thread>
#include <chrono>
#include <iostream>

void foo() {
	std::this_thread::sleep_for(std::chrono::seconds(3));
}

int main() {
	std::thread t;
	t = std::thread(foo);
	for (;;){
		std::cout << t.joinable(<< std::endl;
	}
}
