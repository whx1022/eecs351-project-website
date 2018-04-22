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
    from: '2017-01-01',
    to: '2018-01-01',
    pchart: null,
    range: [],
  },
  watch:{
    sma_n: function(){
      this.drawPrice();
    },
    ewma_a: function(){
      this.drawPrice();
    },
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

    test: function () {
      var price=[];
      var date=[];
      var self=this;
      axios.get('https://www.alphavantage.co/query?function=TIME_SERIES_DAILY_ADJUSTED&&symbol='+self.symbol+'&outputsize=full&apikey=R2JMI1FJK40J7G3G')
        .then(function(res) {
          for(var key in res.data['Time Series (Daily)']){
            if (key>=self.from && key<=self.to){
              price.push(Number(res.data['Time Series (Daily)'][key]['4. close']));
              date.push(new Date(key));
            }
          }
          self.price=price.reverse();
          self.date=date.reverse();
          self.drawPrice();
        })
    },
    drawPrice: function(){
      var data = new google.visualization.DataTable();
      data.addColumn('date', 'Time');
      data.addColumn('number', 'Price');
      data.addColumn('number', 'SMA');
      data.addColumn('number', 'EWMA');
      data.addColumn('number', 'Kalman');
      var sma=this.sma;
      var ewma=this.ewma;
      var kalman=this.kalman;
      for (var i = 0; i < this.price.length; i++) {
          data.addRow([this.date[i],this.price[i],sma[i],ewma[i],kalman[i]])
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

  }
});

