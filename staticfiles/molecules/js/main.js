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

function setMoleculeCookies(name, smiles) {
	const expiry = new Date();
	expiry.setTime(expiry.getTime() + (7 * 24 * 60 * 60 * 1000));
	const expires = `; path=/; expires=${expiry.toUTCString()}`;
	document.cookie = `moleculeSmiles=${smiles}${expires}`;
	document.cookie = `moleculeName=${name}${expires}`;
}

function OpenModel(name, smiles) {
	setMoleculeCookies(name, smiles)

	window.location.href = window.location.origin + '/view_molecule/'
}

document.getElementById("sendNameBtn").addEventListener("click", function () {
	const inputName = document.getElementById("inputName");
	const loader = document.getElementById("nameLoader");
	const loader_container = document.getElementById("loader_container");

	const name = inputName.value.trim();
	if (!name) return;

	loader.style.display = "inline-block";
	loader_container.style.display = "flex";

	axios.post("/send-name/", { name: name }, {
		headers: {
			"X-CSRFToken": getCookie("csrftoken"),
			"Content-Type": "application/json"
		}
	})
	.then(res => {
		const data = res.data;
		if (data.status === "success") {
			setMoleculeCookies(data.name, data.smiles);
			document.getElementById("sendNameBtn").style.display = "none";

			const link = document.getElementById("moleculeLinkName");
			link.href = "/view_molecule/";
			link.style.display = "inline";
		} else {
			Swal.fire({
				title: 'Помилка!',
				text: data.error || 'Щось пішло не так',
				icon: 'error',
				confirmButtonText: 'ОК'
			});
		}
	})
	.catch(error => {
		Swal.fire({
			title: 'Помилка!',
			text: error.response?.data?.error || 'Сталася невідома помилка',
			icon: 'error',
			confirmButtonText: 'ОК'
		});
	})
	.finally(() => {
		loader.style.display = "none";
		loader_container.style.display = "none";
	});
});


document.getElementById("sendFormulaBtn").addEventListener("click", function () {
	const inputFormula = document.getElementById("inputFormula");
	const loader = document.getElementById("formulaLoader");
	const loader_container = document.getElementById("loader_container");

	const formula = inputFormula.value.trim();
	if (!formula) return;

	loader.style.display = "inline-block";
	loader_container.style.display = "flex";

	axios.post("/send-formula/", { formula: formula }, {
		headers: {
			"X-CSRFToken": getCookie("csrftoken"),
			"Content-Type": "application/json"
		}
	})
	.then(res => {
		const data = res.data;
		if (data.status === "success") {
			setMoleculeCookies(data.name, data.smiles);
			document.getElementById("sendFormulaBtn").style.display = "none";

			const link = document.getElementById("moleculeLinkFormula");
			link.href = "/view_molecule/";
			link.style.display = "inline";
		} else {
			Swal.fire({
				title: 'Помилка!',
				text: data.error || 'Щось пішло не так',
				icon: 'error',
				confirmButtonText: 'ОК'
			});
		}
	})
	.catch(error => {
		Swal.fire({
			title: 'Помилка!',
			text: error.response?.data?.error || 'Сталася невідома помилка',
			icon: 'error',
			confirmButtonText: 'ОК'
		});
	})
	.finally(() => {
		loader.style.display = "none";
		loader_container.style.display = "none";
	});
});


var value = 'flex';

document.getElementById("nameOfSubstanceBtn").addEventListener("click", function () {
	const nameOfSubstance = document.getElementById("nameOfSubstance");
	const formulaOfSubstance = document.getElementById("formulaOfSubstance");
	const input_data_container = document.getElementById("input_data_container");


	if (nameOfSubstance.style.display === value) {
		nameOfSubstance.style.display = "none";
		input_data_container.style.display = "none";
	} else {
		nameOfSubstance.style.display = value;
		input_data_container.style.display = value;
		formulaOfSubstance.style.display = "none";
	}
});

document.getElementById("formulaOfSubstanceBtn").addEventListener("click", function () {
	const formulaOfSubstance = document.getElementById("formulaOfSubstance");
	const nameOfSubstance = document.getElementById("nameOfSubstance");
	const input_data_container = document.getElementById("input_data_container");

	if (formulaOfSubstance.style.display === value) {
		formulaOfSubstance.style.display = "none";
		input_data_container.style.display = "none";
	} else {
		formulaOfSubstance.style.display = value;
		input_data_container.style.display = value;
		nameOfSubstance.style.display = "none";
	}
});

document.addEventListener("DOMContentLoaded", function () {
	const container = document.getElementById("molecule-container");
	const emptyMessage = document.getElementById("empty-message");
	const loadMoreButton = document.getElementById("load-more");

	let allMolecules = [];
	let currentIndex = 0;
	const batchSize = 9;

	fetch(window.location.origin + "/api/molecule-history")
		.then(response => response.json())
		.then(data => {
			if (data.length === 0) {
				emptyMessage.style.display = "block";
				return;
			}

			allMolecules = data;
			loadMoreButton.style.display = "inline-block";
			renderNextBatch();
		});

	loadMoreButton.addEventListener("click", renderNextBatch);

	function renderNextBatch() {
		const nextBatch = allMolecules.slice(currentIndex, currentIndex + batchSize);

		nextBatch.forEach(molecule => {
			const elementDiv = document.createElement("div");
			elementDiv.className = "element";
			elementDiv.setAttribute("onclick", `OpenModel(name='${molecule.name}', smiles='${molecule.smiles}')`);

			const nameSpan = document.createElement("span");
			nameSpan.className = "nameOfMol";
			nameSpan.textContent = molecule.name;

			const dateSpan = document.createElement("span");
			dateSpan.className = "date";
			const date = new Date(molecule.created_at);
			const formattedDate = date.toLocaleDateString("uk-UA", {
				year: "2-digit",
				month: "2-digit",
				day: "2-digit"
			}).replace(/\//g, '.');

			dateSpan.textContent = formattedDate;

			elementDiv.appendChild(nameSpan);
			elementDiv.appendChild(dateSpan);
			container.appendChild(elementDiv);
		});

		currentIndex += batchSize;

		if (currentIndex >= allMolecules.length) {
			loadMoreButton.style.display = "none";
		}
	}
});