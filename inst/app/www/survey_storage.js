// survey_storage.js
// Handle survey completion tracking via localStorage

$(document).on('shiny:connected', function() {
  // Check if survey was already completed
  const done = localStorage.getItem('surveyCompleted') === 'true';
  console.log('Survey status on load:', done);
  console.log('Shiny is ready, sending value to server');
  
  // Send status to Shiny using the namespaced input ID
  if (typeof surveyInputId !== 'undefined') {
    Shiny.setInputValue(surveyInputId, done, {priority: 'event'});
  }
});

// Listen for server message to mark survey as complete
Shiny.addCustomMessageHandler('markSurveyComplete', function(msg) {
  console.log('Marking survey complete');
  localStorage.setItem('surveyCompleted', 'true');
  
  // Update Shiny input
  if (typeof surveyInputId !== 'undefined') {
    Shiny.setInputValue(surveyInputId, true, {priority: 'event'});
  }
});
