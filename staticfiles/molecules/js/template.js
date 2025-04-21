document.addEventListener('DOMContentLoaded', function() {
    var burgerMenu = document.getElementById('burger-menu');
    var overlay = document.getElementById('menu');
    var body = document.body;

    burgerMenu.addEventListener('click', function() {
        this.classList.toggle('close');
        overlay.classList.toggle('overlay');

        // Додаємо або прибираємо клас no-scroll для блокування скролінгу
        if (overlay.classList.contains('overlay')) {
            body.classList.add('no-scroll');
        } else {
            body.classList.remove('no-scroll');
        }
    });
});
