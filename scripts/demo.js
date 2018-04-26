google.charts.load('current', {packages: ['corechart', 'line']});

var app;
app = new Vue({
  el: '#demo',
  data: {
    symbol: '',
    price: [],
    sma_n: 10,
    ewma_a: 33,
    date: [],
    from: '',
    to: '',
    actual: [],
    pchart: null,
    rchart: null,
    kchart: null,
    ychart: null,
    range: [],
    p: 4,
    z: 4,
    mse: [0,0,0],
    filter: '',
  },
  watch:{
    sma_n: function(){
      this.drawAll();
    },
    ewma_a: function(){
      this.drawAll();
    },
    p: function(){
      this.drawProny();
    },
    z: function(){
      this.drawProny();
    },
    filter: function(){
      this.drawRegression();
      this.drawKalman();
      this.drawProny();
    }
  },
  computed: {
    sma: function () {
      var sma = [];
      var n=this.sma_n;
      sma[0]=this.price[0]/n;
      for (var i = 1; i < this.price.length; i++) {
        if (i-n>=0){
          sma[i]=sma[i - 1]+ this.price[i] / n-this.price[i - n] / n
        }
        else {
          sma[i]=sma[i - 1]+ this.price[i] / n
        }
      }
      for (var i = 0; i<n;i++){
        sma[i]=NaN;
      }
      return sma;
    },
    ewma: function () {
      var ewma = [];
      var a=this.ewma_a/100;
      ewma[0]=this.price[0];
      for (var i = 1; i < this.price.length; i++) {
        ewma[i]=ewma[i-1]+a*(this.price[i]-ewma[i-1]);
      }
      return ewma;
    },
    kalman: function(){
      var Q=1;
      var R=1;
      var T=1;
      var F=math.matrix([[1,T],[0,1]]);
      var G=math.matrix([[T^2/2],[T]]);
      var H=math.matrix([[1,0]]);
      var size=this.price.length;
      var data=math.zeros(2,size+1);
      var price=math.matrix(this.price);
      data.subset(math.index(0,math.range(1,size+1)),price);
      var delta1=data.subset(math.index(0,math.range(1,size+1)));
      var delta2=data.subset(math.index(0,math.range(0,size+0)));
      var delta=math.subtract(delta1,delta2);
      data.subset(math.index(1,math.range(0,size)),delta);
      var x=math.mean(data,1);
      var xx=math.matrix();
      xx.subset(math.index(0,math.range(0,size+1)),math.multiply(math.ones(size+1),x.subset(math.index(0))));
      xx.subset(math.index(1,math.range(0,size+1)),math.multiply(math.ones(size+1),x.subset(math.index(1))));
      var va=math.subtract(data,xx);
      var v=null;
      var p=math.multiply(va,math.transpose(va));
      var k=null;
      var kalman=[];
      for (var i=0;i<size;i++){
        x=math.multiply(F,x);
        v=math.subtract(data.subset(math.index(0,i+1)),math.multiply(H,x));
        p=math.add(math.multiply(F,p,math.transpose(F)),math.multiply(G,Q,math.transpose(G)));
        k=math.multiply(p,math.transpose(H),math.inv(math.add(math.multiply(H,p,math.transpose(H)),R)));
        x=math.add(x,math.multiply(k,v));
        p=math.multiply(math.subtract(math.eye(2),math.multiply(k,H)),p);
        kalman[i]=x.subset(math.index(0));
      }
      return kalman;
    }
  },
  methods: {
    regression: function(x,y){
      x = x.map(function(d){ return [+d]; });
      y = y.map(function(d){ return [+d]; });
      var m = y.length;
      x = math.concat(math.ones(m,1), x);
      y = math.matrix(y);
      var tr = math.transpose(x);
      var tr_x = math.multiply(tr, x);
      var tr_y = math.multiply(tr, y);
      var theta = math.multiply( math.inv(tr_x), tr_y );
      return function() {
        var args = Array.prototype.slice.call(arguments);
        return math.multiply(math.matrix(math.concat([1], args)), theta)._data[0];
      }
    },
    prony: function(x,length){
      var p=Number(this.p);
      var z=Number(this.z);
      var r=[[0]];
      for (var k=0;k<=p;k++){
        r[k]=[];
        for (var l=0;l<=p;l++){
          r[k][l]=0;
          for (var n=z+1;n<x.length;n++){
            if (n-l>=0 && n-k>=0)
              r[k][l]+=x[n-l]*x[n-k];
          }
        }
      }
      var rm=math.matrix(r);
      var b=math.multiply(-1,math.subset(rm,math.index(math.range(1,p+1),0)));
      var a=math.subset(rm,math.index(math.range(1,p+1),math.range(1,p+1)));
      var den=[1];
      for (k=0;k<p;k++){
        dx=math.subset(a,math.index(math.range(0,p),k),b);
        den[k+1]=math.det(dx) / math.det(a);
      }
      var num=[];
      for (k=0;k<=z;k++){
        num[k]=x[k];
        for (l=1;l<=p;l++){
          if (k-l>=0)
            num[k]+=den[l]*x[k-l];
        }
      }
      var delta=function(x){
        return x==0;
      };
      var y=[];
      for (var i=0;i<x.length+length;i++){
        y[i]=0;
        for (k=0;k<=z;k++){
          y[i]+=delta(i-k)*num[k];
        }
        for (k=1;k<=p;k++){
          if(i-k>=0){
            y[i]-=y[i-k]*den[k];
          }
        }
      }
      return y;
    },
    pkalman: function(X){
      var Q=1;
      var R=1;
      var T=1;
      var F=math.matrix([[1,T],[0,1]]);
      var G=math.matrix([[T^2/2],[T]]);
      var H=math.matrix([[1,0]]);
      var size=X.length;
      var data=math.zeros(2,size+1);
      var price=math.matrix(X);
      data.subset(math.index(0,math.range(1,size+1)),price);
      var delta1=data.subset(math.index(0,math.range(1,size+1)));
      var delta2=data.subset(math.index(0,math.range(0,size+0)));
      var delta=math.subtract(delta1,delta2);
      data.subset(math.index(1,math.range(0,size)),delta);
      var x=math.mean(data,1);
      var xx=math.matrix();
      xx.subset(math.index(0,math.range(0,size+1)),math.multiply(math.ones(size+1),x.subset(math.index(0))));
      xx.subset(math.index(1,math.range(0,size+1)),math.multiply(math.ones(size+1),x.subset(math.index(1))));
      var va=math.subtract(data,xx);
      var v=null;
      var p=math.multiply(va,math.transpose(va));
      var k=null;
      var kalman=[];
      for (var i=0;i<size;i++){
        x=math.multiply(F,x);
        kalman[i]=x.subset(math.index(0));
        v=math.subtract(data.subset(math.index(0,i+1)),math.multiply(H,x));
        p=math.add(math.multiply(F,p,math.transpose(F)),math.multiply(G,Q,math.transpose(G)));
        k=math.multiply(p,math.transpose(H),math.inv(math.add(math.multiply(H,p,math.transpose(H)),R)));
        x=math.add(x,math.multiply(k,v));
        p=math.multiply(math.subtract(math.eye(2),math.multiply(k,H)),p);
      }
      return kalman;
    },
    getPrice: function () {
      if (this.from==='' || this.to==='' || this.symbol===''){
        return;
      }

      var price=[];
      var date=[];
      var actual=[];
      var self=this;
      axios.get('https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&&symbol='+self.symbol+'&outputsize=full&apikey=R2JMI1FJK40J7G3G')
        .then(function(res) {
          for(var key in res.data['Time Series (Daily)']){
            if (key>=self.from && key<=self.to){
              price.push(Number(res.data['Time Series (Daily)'][key]['4. close']));
              date.push(new Date(key));
            }
            else if (key > self.to){
              actual.push(Number(res.data['Time Series (Daily)'][key]['4. close']));
              date.push(new Date(key));
            }
          }
          self.price=price.reverse();
          self.date=date.reverse();
          self.actual=actual.reverse();
          setTimeout(function(){
            self.drawAll();
          },5);
        })
      /*var price=[62.7900000000000,62.4800000000000,62.1900000000000,62.3000000000000,62.7600000000000,62.7300000000000,62.6100000000000,63.0600000000000,62.6200000000000,62.6800000000000,62.6700000000000,62.2400000000000,62.6700000000000,62.7000000000000,63.2000000000000,63.9500000000000,64.1200000000000,65.3900000000000,65.6900000000000,64.8600000000000,64.3550000000000,63.2500000000000,63.5000000000000,63.5000000000000,63.7400000000000,63.5700000000000,63.5200000000000,64.2500000000000,64.2400000000000,64.4100000000000,64.5000000000000,64.7400000000000,64.4700000000000,64.6100000000000,64.3300000000000,64.4200000000000,64.5300000000000,64.5400000000000,64.0800000000000,64.1300000000000,64.6900000000000,63.9900000000000,63.9700000000000,64.1900000000000,64.2600000000000,65.1900000000000,65.1100000000000,65.0100000000000,64.5300000000000,64.5500000000000,64.7500000000000,64.9100000000000,64.9100000000000,65.1900000000000,64.1200000000000,64.9400000000000,65.3600000000000,64.6300000000000,64.9600000000000,65.1200000000000,65.4200000000000,65.6500000000000,65.8100000000000,65.3900000000000,66.3000000000000,65.6000000000000,65.8500000000000,65.6100000000000,65.6000000000000,65.4200000000000,65.2900000000000,65.0400000000000,65.3300000000000,65.6500000000000,65.4600000000000,65.6700000000000];
      var price=[1,4,5,3,2,5,6,3,2,4,5,3,5,4,0.250000000000000,4,6.25000000000000,2.25000000000000,1,6.25000000000000,9,2.25000000000000,1,4,6.25000000000000,2.25000000000000,6.25000000000000,4];
      var date=['03-Jan-2017','04-Jan-2017','05-Jan-2017','06-Jan-2017','09-Jan-2017','10-Jan-2017','11-Jan-2017','12-Jan-2017','13-Jan-2017','17-Jan-2017','18-Jan-2017','19-Jan-2017','20-Jan-2017','23-Jan-2017','24-Jan-2017','25-Jan-2017','26-Jan-2017','27-Jan-2017','30-Jan-2017','31-Jan-2017','01-Feb-2017','02-Feb-2017','03-Feb-2017','06-Feb-2017','07-Feb-2017','08-Feb-2017','09-Feb-2017','10-Feb-2017','13-Feb-2017','14-Feb-2017','15-Feb-2017','16-Feb-2017','17-Feb-2017','21-Feb-2017','22-Feb-2017','23-Feb-2017','24-Feb-2017','27-Feb-2017','28-Feb-2017','01-Mar-2017','02-Mar-2017','03-Mar-2017','06-Mar-2017','07-Mar-2017','08-Mar-2017','09-Mar-2017','10-Mar-2017','13-Mar-2017','14-Mar-2017','15-Mar-2017','16-Mar-2017','17-Mar-2017','20-Mar-2017','21-Mar-2017','22-Mar-2017','23-Mar-2017','24-Mar-2017','27-Mar-2017','28-Mar-2017','29-Mar-2017','30-Mar-2017','31-Mar-2017','03-Apr-2017','04-Apr-2017','05-Apr-2017','06-Apr-2017','07-Apr-2017','10-Apr-2017','11-Apr-2017','12-Apr-2017','13-Apr-2017','17-Apr-2017','18-Apr-2017','19-Apr-2017','20-Apr-2017','21-Apr-2017','24-Apr-2017','25-Apr-2017','26-Apr-2017','27-Apr-2017','28-Apr-2017','01-May-2017','02-May-2017','03-May-2017','04-May-2017','05-May-2017','08-May-2017','09-May-2017','10-May-2017','11-May-2017','12-May-2017','15-May-2017','16-May-2017','17-May-2017','18-May-2017','19-May-2017','22-May-2017','23-May-2017','24-May-2017','25-May-2017','26-May-2017','30-May-2017','31-May-2017','01-Jun-2017','02-Jun-2017','05-Jun-2017','06-Jun-2017','07-Jun-2017','08-Jun-2017','09-Jun-2017','12-Jun-2017','13-Jun-2017','14-Jun-2017','15-Jun-2017','16-Jun-2017','19-Jun-2017','20-Jun-2017','21-Jun-2017','22-Jun-2017','23-Jun-2017','26-Jun-2017','27-Jun-2017','28-Jun-2017','29-Jun-2017','30-Jun-2017','03-Jul-2017','05-Jul-2017','06-Jul-2017','07-Jul-2017','10-Jul-2017','11-Jul-2017','12-Jul-2017','13-Jul-2017','14-Jul-2017','17-Jul-2017','18-Jul-2017','19-Jul-2017','20-Jul-2017','21-Jul-2017','24-Jul-2017','25-Jul-2017','26-Jul-2017','27-Jul-2017','28-Jul-2017','31-Jul-2017','01-Aug-2017','02-Aug-2017','03-Aug-2017','04-Aug-2017','07-Aug-2017','08-Aug-2017','09-Aug-2017','10-Aug-2017','11-Aug-2017','14-Aug-2017','15-Aug-2017','16-Aug-2017','17-Aug-2017','18-Aug-2017','21-Aug-2017','22-Aug-2017','23-Aug-2017','24-Aug-2017','25-Aug-2017','28-Aug-2017','29-Aug-2017','30-Aug-2017','31-Aug-2017','01-Sep-2017','05-Sep-2017','06-Sep-2017','07-Sep-2017','08-Sep-2017','11-Sep-2017','12-Sep-2017','13-Sep-2017','14-Sep-2017','15-Sep-2017','18-Sep-2017','19-Sep-2017','20-Sep-2017','21-Sep-2017','22-Sep-2017','25-Sep-2017','26-Sep-2017','27-Sep-2017','28-Sep-2017','29-Sep-2017','02-Oct-2017','03-Oct-2017','04-Oct-2017','05-Oct-2017','06-Oct-2017','09-Oct-2017','10-Oct-2017','11-Oct-2017','12-Oct-2017','13-Oct-2017','16-Oct-2017','17-Oct-2017','18-Oct-2017','19-Oct-2017','20-Oct-2017','23-Oct-2017','24-Oct-2017','25-Oct-2017','26-Oct-2017','27-Oct-2017','30-Oct-2017','31-Oct-2017','01-Nov-2017','02-Nov-2017','03-Nov-2017','06-Nov-2017','07-Nov-2017','08-Nov-2017','09-Nov-2017','10-Nov-2017','13-Nov-2017','14-Nov-2017','15-Nov-2017','16-Nov-2017','17-Nov-2017','20-Nov-2017','21-Nov-2017','22-Nov-2017','24-Nov-2017','27-Nov-2017','28-Nov-2017','29-Nov-2017','30-Nov-2017','01-Dec-2017','04-Dec-2017','05-Dec-2017','06-Dec-2017','07-Dec-2017','08-Dec-2017','11-Dec-2017','12-Dec-2017','13-Dec-2017','14-Dec-2017','15-Dec-2017','18-Dec-2017','19-Dec-2017','20-Dec-2017','21-Dec-2017','22-Dec-2017','26-Dec-2017','27-Dec-2017','28-Dec-2017','29-Dec-2017'];
      for (var i=0;i<price.length;i++){
        price[i]=Number(price[i]);
      }
      for (var i=0;i<date.length;i++){
        date[i]=new Date(date[i]);
      }
      this.price=price;
      this.date=date;
      this.drawAll();*/
    },
    drawPrice: function(){
      var data = new google.visualization.DataTable();
      data.addColumn('date', 'Time');
      data.addColumn('number', 'Price');
      data.addColumn('number', 'SMA');
      data.addColumn('number', 'EWMA');
      data.addColumn('number', 'Kalman');
      if (this.price.length){
        var sma=this.sma;
        var ewma=this.ewma;
        var kalman=this.kalman;
        for (var i = 0; i < this.price.length; i++) {
          data.addRow([this.date[i],this.price[i],sma[i],ewma[i],kalman[i]])
        }
      }

      var options = {
        hAxis: {
          title: 'Time'
        },
        vAxis: {
          title: 'Price'
        },
        height: 500,
        title:'Raw Stock Price and Filtered Price of ' + this.symbol,
      };
      this.pchart = new google.visualization.LineChart(document.getElementById('price'));
      this.pchart.draw(data, options);
    },
    drawRegression: function(){
      var data = new google.visualization.DataTable();
      data.addColumn('date', 'Time');
      data.addColumn('number', 'Filtered Price');
      data.addColumn('number', 'Actual Price');
      data.addColumn('number', 'Regression');
      data.addColumn('number', 'Prediction');
      if (this.price.length){
        var ts=[];
        for (i=0;i<this.price.length;i++){
          ts[i]=this.date[i].getTime();
        }
        var date=this.date;
        var length=10;
        var filtered=[];
        if (this.filter==='ewma'){
          filtered=this.ewma;
        }
        else if (this.filter==='sma'){
          filtered=this.sma;
          filtered=filtered.slice(this.sma_n,filtered.length);
          ts=ts.slice(this.sma_n,ts.length);
          date=date.slice(this.sma_n,date.length);
        }
        else if (this.filter==='kalman'){
          filtered=this.kalman;
        }
        if (filtered.length){
          var f=this.regression(ts.slice(ts.length-filtered.length,ts.length),filtered);
          var mse=0.;
          for (var i = 0; i < ts.length; i++) {
            data.addRow([date[i],filtered[i],NaN,f(ts[i]),NaN])
          }
          for (; i < ts.length+length; i++) {
            mse+=(this.actual[i-ts.length]-f(date[i].getTime()))*(this.actual[i-ts.length]-f(date[i].getTime()));
            data.addRow([date[i],NaN,this.actual[i-ts.length],NaN,f(date[i].getTime())]);
          }
          console.log(mse);
          mse=mse/length;
          this.mse[0]=mse;
        }
      }
      var options = {
        hAxis: {
          title: 'Time'
        },
        vAxis: {
          title: 'Price'
        },
        height: 500,
      };
      this.rchart = new google.visualization.LineChart(document.getElementById('regression'));
      this.rchart.draw(data, options);
    },
    drawProny: function(){
      var data = new google.visualization.DataTable();
      data.addColumn('date', 'Time');
      data.addColumn('number', 'Filtered Price');
      data.addColumn('number', 'Actual Price');
      data.addColumn('number', 'Prony');
      data.addColumn('number', 'Prediction');
      if (this.price.length){
        var length=10;
        var filtered;
        var date=this.date.slice(0,this.price.length+length);
        if (this.filter==='ewma'){
          filtered=this.ewma;
        }
        else if (this.filter==='sma'){
          filtered=this.sma;
          filtered=filtered.slice(this.sma_n,filtered.length);
          date=date.slice(this.sma_n,date.length);
        }
        else {
          filtered=this.kalman;
        }
        var mse=0.;
        var prony=this.prony(filtered,length);
        for (var i=0; i < date.length-length; i++) {
          data.addRow([date[i],filtered[i],NaN,prony[i],NaN]);
        }
        for (; i < date.length; i++) {
          mse+=Number(this.actual[i-date.length+length]-prony[i])*(this.actual[i-date.length+length]-prony[i]);
          console.log(mse);
          data.addRow([date[i],NaN,this.actual[i-date.length+length],NaN,prony[i]]);
        }
        mse=mse/length;
        this.mse[1]=mse;
      }
      var options = {
        hAxis: {
          title: 'Time'
        },
        vAxis: {
          title: 'Price'
        },
        height: 500,
      };
      this.ychart = new google.visualization.LineChart(document.getElementById('prony'));
      this.ychart.draw(data, options);

    },
    drawKalman: function(){
      var data = new google.visualization.DataTable();
      data.addColumn('date', 'Time');
      data.addColumn('number', 'Actual Price');
      data.addColumn('number', 'Prediction');
      var mse=0.;
      if (this.price.length){
        var length=10;
        var kalman=this.pkalman(this.price.concat(this.actual));
        for (var i = 0; i < this.price.length; i++) {
          data.addRow([this.date[i],this.price[i],kalman[i]]);
        }
        for (; i < this.price.length+length; i++) {
          mse+=(this.actual[i-this.price.length]-kalman[i])*(this.actual[i-this.price.length]-kalman[i]);
          data.addRow([this.date[i],this.actual[i-this.price.length],kalman[i]]);
        }
        mse=mse/length;
        this.mse[2]=mse;
      }
      var options = {
        hAxis: {
          title: 'Time'
        },
        vAxis: {
          title: 'Price'
        },
        height: 500,
      };
      this.kchart = new google.visualization.LineChart(document.getElementById('kalman'));
      this.kchart.draw(data, options);
    },
    drawAll: function(){
      this.drawPrice();
      this.drawRegression();
      this.drawProny();
      this.drawKalman();
    },
    calculateMSE: function(){
      /*var length=10;
      var ts=[];
      for (i=0;i<this.price.length;i++){
        ts[i]=this.date[i].getTime();
      }
      var ewma=this.ewma;
      var kalman=this.ewma;
      var sma=this.sma;
      var fe=this.regression(ts,ewma);
      var sum=0;
      for (var i=0;i<this.actual.length;i++){
        sum+=(this.actual[i]-fe[i+this.price.length])^2;
      }
      console.log(sum);
      this.mse['l']['e']=sum/i;
      var fk=this.regression(ts,kalman);
      sum=0;
      for (var i=0;i<this.actual.length;i++){
        sum+=(this.actual[i]-fk[i+this.price.length])^2;
      }
      console.log(sum);
      this.mse['l']['k']=sum/i;
      var fs=this.regression(ts.slice(ts.length-sma.length,ts.length),sma);
      sum=0;
      for (var i=0;i<this.actual.length;i++){
        sum+=(this.actual[i]-fs[i+this.price.length])^2;
      }
      console.log(sum);
      this.mse['l']['s']=sum/i;
      console.log(this.mse);*/
    }
  },
  mounted: function () {
    var self=this;
    setTimeout(function(){self.drawAll()}, 2000);

  },
});


var from = new mdDateTimePicker.default({
  type: 'date',
  future: moment(),
});
var to = new mdDateTimePicker.default({
  type: 'date',
  future: moment()
});
var fromInput = document.getElementById('from');
fromInput.addEventListener('click', function() {
  from.toggle();
});
from.trigger = document.getElementById('from');
document.getElementById('from').addEventListener('onOk', function() {
  app.from=from.time.format("YYYY-MM-DD");
  to._future= moment.min(moment(),from.time.add(30, 'days'));
  if (moment(document.getElementById('to').value) > to._future){
    app.to=to._future.format("YYYY-MM-DD");
  }
});
var toInput = document.getElementById('to');
toInput.addEventListener('click', function() {
  to.toggle();
});
to.trigger = document.getElementById('to');
document.getElementById('to').addEventListener('onOk', function() {
  app.to=to.time.format("YYYY-MM-DD");
});



document.getElementById('kalmantab').onclick = function(){
  setTimeout(function(){app.drawKalman()}, 5);

  document.getElementById("ffkalman").parentNode.MaterialRadio.check();
  document.getElementById("fsma").parentNode.MaterialRadio.uncheck();
  document.getElementById("fewma").parentNode.MaterialRadio.uncheck();
  document.getElementById("fsma").parentNode.MaterialRadio.disable();
  document.getElementById("fewma").parentNode.MaterialRadio.disable();
};
document.getElementById('pronytab').onclick = function(){
  setTimeout(function(){app.drawProny()}, 5);
  document.getElementById("fsma").parentNode.MaterialRadio.enable();
  document.getElementById("fewma").parentNode.MaterialRadio.enable();
};
document.getElementById('regtab').onclick = function(){
  setTimeout(function(){app.drawRegression()}, 5);
  document.getElementById("fsma").parentNode.MaterialRadio.enable();
  document.getElementById("fewma").parentNode.MaterialRadio.enable();
};
