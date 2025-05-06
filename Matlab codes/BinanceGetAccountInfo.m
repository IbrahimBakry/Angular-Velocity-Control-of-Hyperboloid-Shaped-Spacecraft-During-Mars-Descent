function BinanceGetAccountInfo()

% [key,secret]=key_secret('binance');  
key ='vFXO7ukjUY4e5K2Pk7OvhizyWBFm3aEPCfRfSnZRLdM5kCZLgNOa41E6wcPZZrQL';
secret ='E9MShS7IOwcIuU7HzcxPw5fs6LREuqGCfGQQsd9qA1gioQnKROukNZYsTW1K10HY';
timestamp=binanceServerTime;

queryString =['timestamp=' timestamp]
Signature = crypto(queryString, secret, 'HmacSHA256');
Signature=string(Signature)

url='https://api.binance.com/api/v3/';
url_ext='account';
url=[url url_ext '?' queryString '&signature=' Signature]

postparams=['X-MBX-APIKEY=' key]

header=http_createHeader('Content-Type','application/x-www-form-urlencoded')


[response,status] = urlread2(url,'POST',postparams,header);
verifStatus=status.status

end

function signStr = crypto(str, key, algorithm)
import java.net.*;
import javax.crypto.*;
import javax.crypto.spec.*;
import org.apache.commons.codec.binary.*

keyStr = java.lang.String(key);
key = SecretKeySpec(keyStr.getBytes('UTF-8'), algorithm);
mac = Mac.getInstance(algorithm);
mac.init(key);
toSignStr = java.lang.String(str);
signStr = java.lang.String(Hex.encodeHex( mac.doFinal(  toSignStr.getBytes('UTF-8'))));
end

function serverTime=binanceServerTime(adTime)
if nargin<1
adTime=0; %millisecondes
end
serverTime=urlread2('https://www.binance.com/api/v1/time');
serverTime=JSON.parse(serverTime);
serverTime=num2str(serverTime.serverTime+adTime);
end