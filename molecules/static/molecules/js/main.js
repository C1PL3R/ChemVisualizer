function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== '') {
      const cookies = document.cookie.split(';');
      for (let i = 0; i < cookies.length; i++) {
          const cookie = cookies[i].trim();
          if (cookie.substring(0, name.length + 1) === (name + '=')) {
              cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
              break;
          }
      }
  }
  return cookieValue;
}

function sendMoleculeRequest(smiles) {
  fetch(window.location.origin + "/molecule-smiles/", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      'X-CSRFToken': getCookie('csrftoken'),
    },
    body: JSON.stringify({ name: smiles.trim() }), // Виправлено ключ
  })
    .then((response) => {
      if (!response.ok) {
        throw new Error(`HTTP Error! Status: ${response.status}`);
      }
      return response.json();
    })
    .then((data) => {
      if (data.status === "success") {
        console.log(data.smiles);
        var link = document.getElementById("moleculeLink");

        link.href = `view_molecule/${data.name}/${data.smiles}`;
        link.style.display = "flex";

      } else {
        console.error(`Помилка сервера: ${data.error}`);
      }
    })
    .catch((error) => {
      console.error("Помилка запиту:", error.message);
      alert("Введіть назву речовини коректно!"); // Виведе сповіщення
    });
}

function SendName() {
  let name = document.getElementById("inputName").value.trim();
  console.log(name)

  if (!name) {
    console.error("Помилка: поле введення порожнє!");
    alert("Будь ласка, введіть назву молекули!");
    return;
  }

  sendMoleculeRequest(encodeURIComponent(name)); // Кодуємо спеціальні символи
}
