#include "simulation.hpp"
#include <SFML/Graphics.hpp>

#include <chrono>
#include <thread>


int main() {
    Simulation simulation("/home/ben/lather/config.cfg");
    double time = 0.0;

    sf::RenderWindow window(sf::VideoMode(1000, 1000), "LATHER display");
    auto image = sf::Image();
    sf::Texture texture;
    texture.create(1000, 1000);

    sf::Sprite sprite(texture);

    std::vector<uint8_t> pixels(1000*1000*4);

    sf::Event event;
    while (window.isOpen()) {
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        auto times = std::vector<double>{time};
        simulation.observe_rv(times, 5000.0, 5001.0);
        auto data = simulation.draw_pixmap(time);
        for (int i = 0; i < data.size(); i++) {
            pixels[i] = data[i];
            pixels[i + 1] = data[i];
            pixels[i + 2] = data[i];
            pixels[i + 3] = 255;
        }

        texture.update(pixels.data());
        window.clear();
        window.draw(sprite);
        window.display();

        time += 0.02;
        std::this_thread::sleep_for(std::chrono::milliseconds(15));
    }

    return 0;
}