$(document).ready( function () {

    $(window).trigger("resize");

    $('#dashboard').DataTable({
      searching: false,
      pageLength: 10,
      lengthChange: false
    });

    $('.gene_report').DataTable({
      searching: true,
      pageLength: 10,
      lengthChange: false
    })
    new DataTable('#iscn', {
      searching: false,
      pageLength: 10,
      lengthChange: false
    });
    $('.segment').DataTable( {
      searching: false,
      pageLength: 10,
      lengthChange: false
    });
    $('.bed').DataTable( {
      searching: false,
      pageLength: 10,
      lengthChange: false
    });
    //TODO fix this part
    // $(".container").hide();
    // $('.report_page_top,.report_page_bottom').bootpag({
    //     total: 10,
    //     page: 1
    // }).on("page", function(event, /* page number here */ num){
    //     $(".container").hide();
    //     $("#"+num).show();
    // });
     
  });