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

// function sendMoleculeRequest(smiles) {
//   fetch(window.location.origin + "/from-name-molecule-get-smiles/", {
//     method: "POST",
//     headers: {
//       "Content-Type": "application/json",
//       'X-CSRFToken': getCookie('csrftoken'),
//     },
//     body: JSON.stringify({ name: smiles.trim() }), // Виправлено ключ
//   })
//     .then((response) => {
//       if (!response.ok) {
//         throw new Error(`HTTP Error! Status: ${response.status}`);
//       }
//       return response.json();
//     })
//     .then((data) => {
//       if (data.status === "success") {
//         console.log(data.smiles);
//         var link = document.getElementById("moleculeLinkName");

//         link.href = `view_molecule/${data.name}/${data.smiles}`;
//         link.style.display = "flex";

//       } else {
//         console.error(`Помилка сервера: ${data.error}`);
//       }
//     })
//     .catch((error) => {
//       console.error("Помилка запиту:", error.message);
//       alert("Введіть назву речовини коректно!"); // Виведе сповіщення
//     });
// }

function SendNameAndSmiles() {
  let nameInput = document.getElementById("inputName");
  let smilesInput = document.getElementById("inputSmiles");
  let link = "";

  let name = nameInput.value.trim();
  let smiles = smilesInput.value.trim();

  if (!name && !smiles) {
    console.error("Помилка: поле введення порожнє!");
    alert("Будь ласка, введіть коректно назву молекули aбо smiles-код!");
    return;
  }

  function setCookie(name, smiles, days) {
    let expires = "";

    if (days) {
      const date = new Date();
      date.setTime(date.getTime() + days * 24 * 60 * 60 * 1000); // Додаємо час життя у мілісекундах
      expires = "; expires=" + date.toUTCString();
    }

    if (name && smiles) {
      alert('Впишіть або назву речовини або smiles!')
    } else if (!name && smiles) {
      document.cookie = `moleculeSmiles=${smiles}${expires}; path=/`;
      link = document.getElementById("moleculeLinkSmileCode");
    } else if (name && !smiles) {
      document.cookie = `moleculeName=${name}${expires}; path=/`;
      link = document.getElementById("moleculeLinkName");
    } else {
      alert('Впишіть або назву речовини або smiles!')
    }
  }

  setCookie(name, smiles, 2);
  link.href = window.location.origin + '/view_molecule/';
  link.style.display = "flex";
  nameInput.value = "";
  smilesInput.value = "";
}


function GenerationSmilesCode() {
  let smiles = document.getElementById("inputSmiles").value;

  if (!smiles) {
    alert("Будь ласка, введіть smiles-код молекули!");
    return;
  }

  fetch(window.location.origin + "/check-smiles-code/", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      'X-CSRFToken': getCookie('csrftoken'),
    },
    body: JSON.stringify({ smiles: smiles }), // Виправлено ключ
  })
    .then((response) => {
      if (!response.ok) {
        throw new Error(`HTTP Error! Status: ${response.status}`);
      }
      return response.json();
    })
    .then((data) => {
      if (data.status === "success") {
        var link = document.getElementById("moleculeLinkSmileCode");

        link.href = `view_molecule-smiles/${smiles}`;
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