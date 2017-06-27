#include "simulation.hpp"
#include <SFML/Graphics.hpp>
#include <chrono>
#include <thread>


int main() {
    Simulation simulation("/home/ben/lather/config.cfg");
    double time = 0.0;

    sf::RenderWindow window(sf::VideoMode(1000, 1000), "LATHER display", sf::Style::Fullscreen);
    window.setVerticalSyncEnabled(true);
    auto image = sf::Image();
    sf::Texture texture;
    texture.create(1000, 1000);
    sf::Sprite sprite(texture);
    sprite.move((1920-1000)/2., (1080-1000)/2.);

    sf::Event event;

    while (window.isOpen()) {

        auto start = std::chrono::high_resolution_clock::now();

        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape) {
                window.close();
            }
        }

        auto pixels = simulation.draw_rgba(time);
        texture.update(pixels.data());

        window.clear();
        window.draw(sprite);
        window.display();

        time += 0.02;

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = end-start;
        std::chrono::milliseconds framerate(33);
        std::this_thread::sleep_for(framerate-duration);
    }

    return 0;
}